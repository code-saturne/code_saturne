# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2011 EDF S.A.
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
Actions Handler
===============

Creates menu, actions, and separators for the SALOME Desktop.
"""

#-------------------------------------------------------------------------------
# Standard modules
#-------------------------------------------------------------------------------

import os
import re
import logging

#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------

from PyQt4.QtGui import QApplication, QCursor, QDialog, QMessageBox
from PyQt4.QtCore import Qt, QObject, SIGNAL, QFileInfo, QString

#-------------------------------------------------------------------------------
# Salome modules
#-------------------------------------------------------------------------------

import VISU
import visu_gui
import SALOMEDS #Si on veut changer de couleur...
import salome
import SMESH
#import smeshDC # smeshDC allows to use smesh with directly accesible methods - for instance GetName()

#-------------------------------------------------------------------------------
# Application modules
#-------------------------------------------------------------------------------

import CFDSTUDYGUI_DialogCollector
import CFDSTUDYGUI_DataModel
import CFDSTUDYGUI_Commons
import CFDSTUDYGUI_CommandMgr
from CFDSTUDYGUI_Agents import *
from CFDSTUDYGUI_Commons import CFD_Code, BinCode, Trace, CFD_Saturne, CFD_Neptune, sgPyQt, sg
import CFDSTUDYGUI_SolverGUI

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("CFDSTUDYGUI_ActionsHandler")
#log.setLevel(logging.DEBUG)
log.setLevel(logging.NOTSET)

#-------------------------------------------------------------------------------
# Global definitions
#-------------------------------------------------------------------------------

#global actions
#CFDSTUDYMenu                 = 0
SetStudyAction                = 1
AddCaseAction                 = 2
RunCaseAction                 = 3
LaunchGUIAction               = 4
OpenXMLCFDGUIAction           = 5
UpdateObjBrowserAction        = 6
InfoCFDSTUDYAction            = 7

#common actions
RemoveAction                  = 20
ViewAction                    = 21
EditAction                    = 22
MoveToDRAFTAction             = 23
CopyInDATAAction              = 24
CopyInSRCAction               = 25
CopyCaseFileAction            = 26

#link actions
CreateLinkAction              = 30
DefineLinkAction              = 31

#export/convert actions
ExportInPostProAction         = 40
ExportInSMESHAction           = 41
ConvertInMEDAction            = 42
ECSConvertAction              = 43

#other actions
CheckCompilationAction        = 50
RunScriptAction               = 51

#Display Actions with VISU
DisplayMESHAction              = 60
DisplayGroupMESHAction         = 61
DisplayOnlyGroupMESHAction     = 62
HideGroupMESHAction            = 63

DisplayTypeMenu                = 70
DisplayTypePOINT               = 71
DisplayTypeWIREFRAME           = 72
DisplayTypeSHADED              = 74
DisplayTypeINSIDEFRAME         = 75
DisplayTypeSURFACEFRAME        = 76
DisplayTypeFEATURE_EDGES       = 77
DisplayTypeSHRINK              = 78

#=====SOLVER ACTIONS
#Common Actions
SolverFileMenu                 = 100
SolverSaveDataFileAction       = 101
SolverSaveAsDataFileAction     = 102

SolverToolsMenu                = 110
SolverOpenShellAction          = 111
SolverDisplayCurrentCaseAction = 112

SolverHelpMenu                 = 130
SolverHelpAboutAction          = 131

#Saturne actions
SaturneReloadModulesAction      = 201
SaturneReloadPageAction         = 202
#Help menu
SaturneHelpLicenseAction        = 251
SaturneHelpUserManualMenu       = 260
SaturneHelpCodeSaturneAction    = 261
SaturneHelpSolutionDomainAction = 262
SaturneHelpCS_KernelAction      = 263
SaturneHelpCS_InfosAction       = 264

#Neptune actions
NeptuneWinMenu                  = 301
NeptuneWinBrowserAction         = 302
NeptuneWinIdentityAction        = 303

StopSolverAction                = 400
ShowSolverProcessAction         = 401

# ObjectTR is a convenient object for traduction purpose

ObjectTR = QObject()

#-------------------------------------------------------------------------------
# Classes definition
#-------------------------------------------------------------------------------

class ActionError(Exception):
    """
    New exception definition.
    """
    def __init__(self, value):
        """
        Constructor.
        """
        self.value = value

    def __str__(self):
        """
        String representation of the attribute I{self.value}.
        """
        return repr(self.value)


class CFDSTUDYGUI_ActionsHandler(QObject):
    def __init__(self):
        """
        Constructor.
        """
        log.debug("__init__")
        QObject.__init__(self, None)

        #intialise all dialogs
        self.DialogCollector = CFDSTUDYGUI_DialogCollector.CFDSTUDYGUI_DialogCollector()

        self._ActionMap = {}
        self._CommonActionIdMap = {}
        self._SolverActionIdMap = {}
        self._SaturneActionIdMap = {}
        self._NeptuneActionIdMap = {}

        self._SalomeSelection = sgPyQt.getSelection()

        self._CommandMgr = CFDSTUDYGUI_CommandMgr.CFDSTUDYGUI_CommandMgr()

        self._SolverGUI = CFDSTUDYGUI_SolverGUI.CFDSTUDYGUI_SolverGUI()

        self._DskAgent = Desktop_Agent()

        self.myVisu = visu_gui.myVisu
        self.myViewManager = self.myVisu.GetViewManager()
        #self.myView = myViewManager.Create3DView()
        #self.myView = myViewManager.GetCurrentView()
        log.debug("__init__ myVisu = %s" % self.myVisu)
        log.debug("__init__ myViewManager = %s" % self.myViewManager)
        #log.debug("__init__ myView = %s" % self.myView)


    def createActions(self):
        """
        Creates menu, actions, and separators.
        """
        menu_id = sgPyQt.createMenu(ObjectTR.tr("CFDSTUDY_MENU"),\
                                     -1,\
                                     -1,\
                                     10)
        tool_id = sgPyQt.createTool(ObjectTR.tr("CFDSTUDY_TOOL_BAR"))

        action = sgPyQt.createAction(-1,\
                                      ObjectTR.tr("SET_CFDSTUDY_STUDY_TEXT"),\
                                      ObjectTR.tr("SET_CFDSTUDY_STUDY_TIP"),\
                                      ObjectTR.tr("SET_CFDSTUDY_STUDY_SB"),\
                                      ObjectTR.tr("SET_CFDSTUDY_STUDY_ICON"))
        sgPyQt.createMenu(action, menu_id)
        sgPyQt.createTool(action, tool_id)

        action_id = sgPyQt.actionId(action)
        self._ActionMap[action_id] = action
        self._CommonActionIdMap[SetStudyAction] = action_id
        self.connect(action, SIGNAL("activated()"), self.slotStudyLocation)

        action = sgPyQt.createSeparator()
        sgPyQt.createMenu(action, menu_id)
        sgPyQt.createTool(action, tool_id)

        action = sgPyQt.createAction(-1,\
                                      ObjectTR.tr("ADD_CFDSTUDY_CASE_TEXT"),\
                                      ObjectTR.tr("ADD_CFDSTUDY_CASE_TIP"),\
                                      ObjectTR.tr("ADD_CFDSTUDY_CASE_SB"),\
                                      ObjectTR.tr("ADD_CFDSTUDY_CASE_ICON"))
        sgPyQt.createMenu(action, menu_id)
        sgPyQt.createTool(action, tool_id)
        action_id = sgPyQt.actionId(action)
        self._ActionMap[action_id] = action
        self._CommonActionIdMap[AddCaseAction] = action_id
        self.connect(action, SIGNAL("activated()"), self.slotAddCase)

        action = sgPyQt.createSeparator()
        sgPyQt.createMenu(action, menu_id)
        sgPyQt.createTool(action, tool_id)

        action = sgPyQt.createAction(-1,\
                                      ObjectTR.tr("LAUNCH_CFDSTUDY_GUI_TEXT"),\
                                      ObjectTR.tr("LAUNCH_CFDSTUDY_GUI_TIP"),\
                                      ObjectTR.tr("LAUNCH_CFDSTUDY_GUI_SB"),\
                                      ObjectTR.tr("LAUNCH_CFDSTUDY_GUI_ICON"))
        sgPyQt.createMenu(action, menu_id)
        sgPyQt.createTool(action, tool_id)
        action_id = sgPyQt.actionId(action)
        self._ActionMap[action_id] = action
        self._CommonActionIdMap[LaunchGUIAction] = action_id
        self.connect(action, SIGNAL("activated()"), self.slotRunGUI)

        action = sgPyQt.createAction(-1,\
                                      ObjectTR.tr("Open CFD GUI"),\
                                      ObjectTR.tr("LAUNCH_CFDSTUDY_GUI_TIP"),\
                                      ObjectTR.tr("LAUNCH_CFDSTUDY_GUI_SB"),\
                                      ObjectTR.tr("LAUNCH_CFDSTUDY_GUI_ICON"))
        sgPyQt.createMenu(action, menu_id)
        sgPyQt.createTool(action, tool_id)
        action_id = sgPyQt.actionId(action)
        self._ActionMap[action_id] = action
        self._CommonActionIdMap[OpenXMLCFDGUIAction] = action_id
        self.connect(action, SIGNAL("activated()"), self.slotOpenCFD_GUI)

        action = sgPyQt.createAction(-1,\
                                      ObjectTR.tr("LAUNCH_CFDSTUDY_CASE_TEXT"),\
                                      ObjectTR.tr("LAUNCH_CFDSTUDY_CASE_TIP"),\
                                      ObjectTR.tr("LAUNCH_CFDSTUDY_CASE_SB"),\
                                      ObjectTR.tr("LAUNCH_CFDSTUDY_CASE_ICON"))
        sgPyQt.createMenu(action, menu_id)
        sgPyQt.createTool(action, tool_id)
        action_id = sgPyQt.actionId(action)
        self._ActionMap[action_id] = action
        self._CommonActionIdMap[RunCaseAction] = action_id
        self.connect(action, SIGNAL("activated()"), self.slotRunCase)

        action = sgPyQt.createSeparator()
        sgPyQt.createMenu(action, menu_id)
        sgPyQt.createTool(action, tool_id)

        action = sgPyQt.createAction(-1,\
                                      ObjectTR.tr("UPDATE_CFDSTUDY_OBJBROWSER_TEXT"),\
                                      ObjectTR.tr("UPDATE_CFDSTUDY_OBJBROWSER_TIP"),\
                                      ObjectTR.tr("UPDATE_CFDSTUDY_OBJBROWSER_SB"),\
                                      ObjectTR.tr("UPDATE_CFDSTUDY_OBJBROWSER_ICON"))
        sgPyQt.createMenu(action, menu_id)
        sgPyQt.createTool(action, tool_id)
        action_id = sgPyQt.actionId(action)
        self._ActionMap[action_id] = action
        self._CommonActionIdMap[UpdateObjBrowserAction] = action_id
        self.connect(action, SIGNAL("activated()"), self.slotUpdateObjectBrowser)

        action = sgPyQt.createSeparator()
        sgPyQt.createMenu(action, menu_id)

        action = sgPyQt.createAction(-1,\
                                      ObjectTR.tr("INFO_CFDSTUDY_TEXT"),\
                                      ObjectTR.tr("INFO_CFDSTUDY_TIP"),\
                                      ObjectTR.tr("INFO_CFDSTUDY_SB"),\
                                      ObjectTR.tr("INFO_CFDSTUDY_ICON"))
        sgPyQt.createMenu(action, menu_id)
        action_id = sgPyQt.actionId(action)
        self._ActionMap[action_id] = action
        self._CommonActionIdMap[InfoCFDSTUDYAction] = action_id
        self.connect(action, SIGNAL("activated()"), self.slotInfo)


        # common actions
        action = sgPyQt.createAction(-1,\
                                      ObjectTR.tr("REMOVE_ACTION_TEXT"),\
                                      ObjectTR.tr("REMOVE_ACTION_TIP"),\
                                      ObjectTR.tr("REMOVE_ACTION_SB"),\
                                      ObjectTR.tr("REMOVE_ACTION_ICON"))
        action_id = sgPyQt.actionId(action)
        self._ActionMap[action_id] = action
        self._CommonActionIdMap[RemoveAction] = action_id
        self.connect(action, SIGNAL("activated()"), self.slotRemoveAction)

        action = sgPyQt.createAction(-1,\
                                      ObjectTR.tr("VIEW_ACTION_TEXT"),\
                                      ObjectTR.tr("VIEW_ACTION_TIP"),\
                                      ObjectTR.tr("VIEW_ACTION_SB"),\
                                      ObjectTR.tr("VIEW_ACTION_ICON"))
        self.connect(action, SIGNAL("activated()"), self.slotViewAction)
        action_id = sgPyQt.actionId(action)
        self._ActionMap[action_id] = action
        self._CommonActionIdMap[ViewAction] = action_id

        action = sgPyQt.createAction(-1,\
                                      ObjectTR.tr("EDIT_ACTION_TEXT"),\
                                      ObjectTR.tr("EDIT_ACTION_TIP"),\
                                      ObjectTR.tr("EDIT_ACTION_SB"),\
                                      ObjectTR.tr("EDIT_ACTION_ICON"))
        self.connect(action, SIGNAL("activated()"), self.slotEditAction)
        action_id = sgPyQt.actionId(action)
        self._ActionMap[action_id] = action
        self._CommonActionIdMap[EditAction] = action_id

        action = sgPyQt.createAction(-1,\
                                      ObjectTR.tr("MOVE_TO_DRAFT_ACTION_TEXT"),\
                                      ObjectTR.tr("MOVE_TO_DRAFT_ACTION_TIP"),\
                                      ObjectTR.tr("MOVE_TO_DRAFT_ACTION_SB"),\
                                      ObjectTR.tr("MOVE_ACTION_ICON"))
        self.connect(action, SIGNAL("activated()"), self.slotMoveToDRAFT)
        action_id = sgPyQt.actionId(action)
        self._ActionMap[action_id] = action
        self._CommonActionIdMap[MoveToDRAFTAction] = action_id

        action = sgPyQt.createAction(-1,\
                                      ObjectTR.tr("COPY_IN_DATA_ACTION_TEXT"),\
                                      ObjectTR.tr("COPY_IN_DATA_ACTION_TIP"),\
                                      ObjectTR.tr("COPY_IN_DATA_ACTION_SB"),\
                                      ObjectTR.tr("COPY_ACTION_ICON"))
        self.connect(action, SIGNAL("activated()"), self.slotCopyInDATA)
        action_id = sgPyQt.actionId(action)
        self._ActionMap[action_id] = action
        self._CommonActionIdMap[CopyInDATAAction] = action_id

        action = sgPyQt.createAction(-1,\
                                     ObjectTR.tr("COPY_IN_SRC_ACTION_TEXT"),\
                                     ObjectTR.tr("COPY_IN_SRC_ACTION_TIP"),\
                                     ObjectTR.tr("COPY_IN_SRC_ACTION_SB"),\
                                     ObjectTR.tr("COPY_ACTION_ICON"))
        self.connect(action, SIGNAL("activated()"), self.slotCopyInSRC)
        action_id = sgPyQt.actionId(action)
        self._ActionMap[action_id] = action
        self._CommonActionIdMap[CopyInSRCAction] = action_id

        action = sgPyQt.createAction(-1,\
                                      ObjectTR.tr("COPY_CASE_FILE_ACTION_TEXT"),\
                                      ObjectTR.tr("COPY_CASE_FILE_ACTION_TIP"),\
                                      ObjectTR.tr("COPY_CASE_FILE_ACTION_SB"),\
                                      ObjectTR.tr("COPY_ACTION_ICON"))
        self.connect(action, SIGNAL("activated()"), self.slotCopyCaseFile)
        action_id = sgPyQt.actionId(action)
        self._ActionMap[action_id] = action
        self._CommonActionIdMap[CopyCaseFileAction] = action_id

        #export/convert actions
        action = sgPyQt.createAction(-1,\
                                      ObjectTR.tr("EXPORT_IN_POSTPRO_ACTION_TEXT"),\
                                      ObjectTR.tr("EXPORT_IN_POSTPRO_ACTION_TIP"),\
                                      ObjectTR.tr("EXPORT_IN_POSTPRO_ACTION_SB"),\
                                      ObjectTR.tr("EXPORT_IN_POSTPRO_ACTION_ICON"))
        self.connect(action, SIGNAL("activated()"), self.slotExportInPostPro)
        action_id = sgPyQt.actionId(action)
        self._ActionMap[action_id] = action
        self._CommonActionIdMap[ExportInPostProAction] = action_id

        action = sgPyQt.createAction(-1,\
                                      ObjectTR.tr("EXPORT_IN_SMESH_ACTION_TEXT"),\
                                      ObjectTR.tr("EXPORT_IN_SMESH_ACTION_TIP"),\
                                      ObjectTR.tr("EXPORT_IN_SMESH_ACTION_SB"),\
                                      ObjectTR.tr("EXPORT_IN_SMESH_ACTION_ICON"))
        self.connect(action, SIGNAL("activated()"), self.slotExportInSMESH)
        action_id = sgPyQt.actionId(action)
        self._ActionMap[action_id] = action
        self._CommonActionIdMap[ExportInSMESHAction] = action_id

# popup added to vizualise the mesh. It is not necessary to switch into SMESH Component

        action = sgPyQt.createAction(-1,\
                                      "Display mesh",\
                                      "Display mesh",\
                                      "Display mesh",\
                                      ObjectTR.tr("MESH_OBJ_ICON"))
        self.connect(action, SIGNAL("activated()"), self.slotDisplayMESH)
        action_id = sgPyQt.actionId(action)
        self._ActionMap[action_id] = action
        self._CommonActionIdMap[DisplayMESHAction] = action_id

# popup added to vizualise the mesh groups. It is not necessary to switch into SMESH Component

        action = sgPyQt.createAction(-1,\
                                      "Display",\
                                      "Display",\
                                      "Display",\
                                      ObjectTR.tr("MESH_TREE_OBJ_ICON"))
        self.connect(action, SIGNAL("activated()"), self.slotDisplayMESHGroups)
        action_id = sgPyQt.actionId(action)
        self._ActionMap[action_id] = action
        self._CommonActionIdMap[DisplayGroupMESHAction] = action_id

        action = sgPyQt.createAction(-1,\
                                      "Display only",\
                                      "Display only",\
                                      "Display only",\
                                      ObjectTR.tr("MESH_TREE_OBJ_ICON"))
        self.connect(action, SIGNAL("activated()"), self.slotDisplayOnlyMESHGroups)
        action_id = sgPyQt.actionId(action)
        self._ActionMap[action_id] = action
        self._CommonActionIdMap[DisplayOnlyGroupMESHAction] = action_id

        action = sgPyQt.createAction(-1,\
                                      "Hide",\
                                      "Hide",\
                                      "Hide",\
                                      ObjectTR.tr("MESH_TREE_OBJ_ICON"))
        self.connect(action, SIGNAL("activated()"), self.slotHideMESHGroups)
        action_id = sgPyQt.actionId(action)
        self._ActionMap[action_id] = action
        self._CommonActionIdMap[HideGroupMESHAction] = action_id

        action = sgPyQt.createAction(-1,\
                                      ObjectTR.tr("CONVERT_IN_MED_ACTION_TEXT"),\
                                      ObjectTR.tr("CONVERT_IN_MED_ACTION_TIP"),\
                                      ObjectTR.tr("CONVERT_IN_MED_ACTION_SB"),\
                                      ObjectTR.tr("CONVERT_IN_MED_ACTION_ICON"))
        action_id = sgPyQt.actionId(action)
        self._ActionMap[action_id] = action
        self._CommonActionIdMap[ConvertInMEDAction] = action_id

        action = sgPyQt.createAction(-1,\
                                      ObjectTR.tr("ECS_CONVERT_ACTION_TEXT"),\
                                      ObjectTR.tr("ECS_CONVERT_ACTION_TIP"),\
                                      ObjectTR.tr("ECS_CONVERT_ACTION_SB"),\
                                      ObjectTR.tr("ECS_CONVERT_ACTION_ICON"))
        self.connect(action, SIGNAL("activated()"), self.slotMeshConvertToMed)
        action_id = sgPyQt.actionId(action)
        self._ActionMap[action_id] = action
        self._CommonActionIdMap[ECSConvertAction] = action_id

        #link actions
        action = sgPyQt.createAction(-1,\
                                      ObjectTR.tr("CREATE_LINK_ACTION_TEXT"),\
                                      ObjectTR.tr("CREATE_LINK_ACTION_TIP"),\
                                      ObjectTR.tr("CREATE_LINK_ACTION_SB"),\
                                      ObjectTR.tr("CREATE_LINK_ACTION_ICON"))
        self.connect(action, SIGNAL("activated()"), self.slotCreateLink)
        action_id = sgPyQt.actionId(action)
        self._ActionMap[action_id] = action
        self._CommonActionIdMap[CreateLinkAction] = action_id

        action = sgPyQt.createAction(-1,\
                                      ObjectTR.tr("DEFINE_LINK_ACTION_TEXT"),\
                                      ObjectTR.tr("DEFINE_LINK_ACTION_TIP"),\
                                      ObjectTR.tr("DEFINE_LINK_ACTION_SB"),\
                                      ObjectTR.tr("DEFINE_LINK_ACTION_ICON"))
        self.connect(action, SIGNAL("activated()"), self.slotDefineLink)
        action_id = sgPyQt.actionId(action)
        self._ActionMap[action_id] = action
        self._CommonActionIdMap[DefineLinkAction] = action_id

        #other actions
        action = sgPyQt.createAction(-1,\
                                      ObjectTR.tr("CHECK_COMPILATION_ACTION_TEXT"),\
                                      ObjectTR.tr("CHECK_COMPILATION_ACTION_TIP"),\
                                      ObjectTR.tr("CHECK_COMPILATION_ACTION_SB"),\
                                      ObjectTR.tr("CHECK_COMPILATION_ACTION_ICON"))
        self.connect(action, SIGNAL("activated()"), self.slotCheckCompilation)
        action_id = sgPyQt.actionId(action)
        self._ActionMap[action_id] = action
        self._CommonActionIdMap[CheckCompilationAction] = action_id

        action = sgPyQt.createAction(-1,\
                                      ObjectTR.tr("RUN_SCRIPT_ACTION_TEXT"),\
                                      ObjectTR.tr("RUN_SCRIPT_ACTION_TIP"),\
                                      ObjectTR.tr("RUN_SCRIPT_ACTION_SB"),\
                                      ObjectTR.tr("RUN_SCRIPT_ACTION_ICON"))
        self.connect(action, SIGNAL("activated()"), self.slotRunScript)
        action_id = sgPyQt.actionId(action)
        self._ActionMap[action_id] = action
        self._CommonActionIdMap[RunScriptAction] = action_id

        # Solver actions

        # File menu
        #Info: Find the menu File into the Main Menu Bar of Salome
        fileId = sgPyQt.createMenu( ObjectTR.tr("MEN_DESK_FILE"), -1,-1)
        #Info: create my menu into  menu File at position 7
        action_id = sgPyQt.createMenu(ObjectTR.tr("SOLVER_FILE_MENU_TEXT"),fileId,-1,7,1)

        #Info: Warning: a Separator is a QMenu item (a trait)
        #Info: create a separator after my menu in position 8
        action = sgPyQt.createSeparator()
        sgPyQt.createMenu(action, fileId,-1,8,1)

        self._SolverActionIdMap[SolverFileMenu] = action_id
        action = sgPyQt.createAction(SolverSaveDataFileAction,\
                                      ObjectTR.tr("SOLVER_SAVE_ACTION_TEXT"),\
                                      ObjectTR.tr("SOLVER_SAVE_ACTION_TIP"),\
                                      ObjectTR.tr("SOLVER_SAVE_ACTION_SB"),\
                                      ObjectTR.tr("SOLVER_SAVE_ACTION_ICON"),
                                      Qt.SHIFT+Qt.CTRL+Qt.Key_S)
        sgPyQt.createTool(action, tool_id)
        sgPyQt.createMenu(action, self._SolverActionIdMap[SolverFileMenu], 100)
        self.connect(action, SIGNAL("activated()"), self.slotSaveDataFile)
        action_id = sgPyQt.actionId(action)
        self._ActionMap[action_id] = action
        self._SolverActionIdMap[SolverSaveDataFileAction] = action_id
        action = sgPyQt.createAction(-1,\
                                      ObjectTR.tr("SOLVER_SAVEAS_ACTION_TEXT"),\
                                      ObjectTR.tr("SOLVER_SAVEAS_ACTION_TIP"),\
                                      ObjectTR.tr("SOLVER_SAVEAS_ACTION_SB"),\
                                      ObjectTR.tr("SOLVER_SAVEAS_ACTION_ICON"),
                                      Qt.SHIFT+Qt.CTRL+Qt.Key_A)
        sgPyQt.createMenu(action, self._SolverActionIdMap[SolverFileMenu], 100)
        self.connect(action, SIGNAL("activated()"), self.slotSaveAsDataFile)
        action_id = sgPyQt.actionId(action)
        self._ActionMap[action_id] = action
        self._SolverActionIdMap[SolverSaveAsDataFileAction] = action_id
        action = sgPyQt.createSeparator()
        sgPyQt.createMenu(action, 1, 0, 2)
        #Tools Menu
        action = sgPyQt.createSeparator()
        sgPyQt.createMenu(action, menu_id, 0, -1)

        action_id = sgPyQt.createMenu(ObjectTR.tr("SOLVER_TOOLS_MENU_TEXT"), menu_id)
        self._SolverActionIdMap[SolverToolsMenu] = action_id
        action = sgPyQt.createAction(-1,\
                                      ObjectTR.tr("SOLVER_OPENSHELL_ACTION_TEXT"),\
                                      ObjectTR.tr("SOLVER_OPENSHELL_ACTION_TIP"),\
                                      ObjectTR.tr("SOLVER_OPENSHELL_ACTION_SB"))
        action_id = sgPyQt.actionId(action)
        self._ActionMap[action_id] = action
        self._SolverActionIdMap[SolverOpenShellAction] = action_id
        self.connect(action, SIGNAL("activated()"), self.slotOpenShell)

        sgPyQt.createMenu(action, self._SolverActionIdMap[SolverToolsMenu])

        action = sgPyQt.createSeparator()
        sgPyQt.createMenu(action, SolverToolsMenu, 0, -1)

        action = sgPyQt.createAction(-1,\
                                      ObjectTR.tr("SOLVER_DISPLAYCASE_ACTION_TEXT"),\
                                      ObjectTR.tr("SOLVER_DISPLAYCASE_ACTION_TIP"),\
                                      ObjectTR.tr("SOLVER_DISPLAYCASE_ACTION_SB"))
        action_id = sgPyQt.actionId(action)
        self._ActionMap[action_id] = action
        self._SolverActionIdMap[SolverDisplayCurrentCaseAction] = action_id
        self.connect(action, SIGNAL("activated()"), self.slotDisplayCurrentCase)

        sgPyQt.createMenu(action, self._SolverActionIdMap[SolverToolsMenu])
        action = sgPyQt.createSeparator()
        sgPyQt.createMenu(action, SolverToolsMenu, 0, -1)
        #for auto hide last separator in tools menu
        self._SaturneActionIdMap[0] = action_id

        # Help menu: insert a Solver Menu Help to the Main Menu Help of Salome

        helpId = sgPyQt.createMenu( ObjectTR.tr("MEN_DESK_HELP"), -1, -1)
        #Info: Separator created at the end of the Menu Help (when we did not indicate a number)

        action = sgPyQt.createSeparator()
        sgPyQt.createMenu(action, helpId)
        #global SolverHelpMenu
        #Info: Solver Help Menu created at the end of the Menu Help of Salome(when we did not indicate a number)
        action_id = sgPyQt.createMenu(ObjectTR.tr("SOLVER_HELP_MENU_TEXT"),helpId)
        self._SolverActionIdMap[SolverHelpMenu] = action_id
        action = sgPyQt.createAction(-1,\
                                     ObjectTR.tr("SOLVER_HELPABOUT_ACTION_TEXT"),\
                                     ObjectTR.tr("SOLVER_HELPABOUT_ACTION_TIP"),\
                                     ObjectTR.tr("SOLVER_HELPABOUT_ACTION_SB"))
        action_id = sgPyQt.actionId(action)
        self._ActionMap[action_id] = action
        self._SolverActionIdMap[SolverHelpAboutAction] = action_id
        sgPyQt.createMenu(action, self._SolverActionIdMap[SolverHelpMenu])
        self.connect(action, SIGNAL("activated()"), self.slotHelpAbout)
        self._ActionMap[action_id].setVisible(True)

        # Saturne actions
        # Tools menu
        action = sgPyQt.createAction(-1,\
                                      ObjectTR.tr("SATURNE_RELOADMODULES_ACTION_TEXT"),\
                                      ObjectTR.tr("SATURNE_RELOADMODULES_ACTION_TIP"),\
                                      ObjectTR.tr("SATURNE_RELOADMODULES_ACTION_SB"))
        sgPyQt.createMenu(action, SolverToolsMenu)
        action_id = sgPyQt.actionId(action)
        self._ActionMap[action_id] = action
        self._SaturneActionIdMap[SaturneReloadModulesAction] = action_id
        self.connect(action, SIGNAL("activated()"), self.slotSaturneReloadModule)

        action = sgPyQt.createAction(-1,\
                                      ObjectTR.tr("SATURNE_RELOADPAGE_ACTION_TEXT"),\
                                      ObjectTR.tr("SATURNE_RELOADPAGE_ACTION_TIP"),\
                                      ObjectTR.tr("SATURNE_RELOADPAGE_ACTION_SB"))
        sgPyQt.createMenu(action, SolverToolsMenu)
        action_id = sgPyQt.actionId(action)
        self._ActionMap[action_id] = action
        self._SaturneActionIdMap[SaturneReloadPageAction] = action_id
        self.connect(action, SIGNAL("activated()"), self.slotSaturneReloadPage)

        # Help menu
        action = sgPyQt.createAction(SaturneHelpLicenseAction,\
                                      ObjectTR.tr("SATURNE_HELPLICENSE_ACTION_TEXT"),\
                                      ObjectTR.tr("SATURNE_HELPLICENSE_ACTION_TIP"),\
                                      ObjectTR.tr("SATURNE_HELPLICENSE_ACTION_SB"))
        sgPyQt.createMenu(action, self._SolverActionIdMap[SolverHelpMenu])
        action_id = sgPyQt.actionId(action)
        self._ActionMap[action_id] = action
        self._SaturneActionIdMap[SaturneHelpLicenseAction] = action_id
        self.connect(action, SIGNAL("activated()"), self.slotSaturneHelpLicense)

        action_id = sgPyQt.createMenu(ObjectTR.tr("SATURNE_USERMANUAL_MENU_TEXT"),\
                                      self._SolverActionIdMap[SolverHelpMenu])
        self._SaturneActionIdMap[SaturneHelpUserManualMenu] = action_id

        action = sgPyQt.createAction(SaturneHelpCodeSaturneAction,\
                                      ObjectTR.tr("SATURNE_HELP_CS_ACTION_TEXT"),\
                                      ObjectTR.tr("SATURNE_HELP_CS_ACTION_TIP"),\
                                      ObjectTR.tr("SATURNE_HELP_CS_ACTION_SB"))
        sgPyQt.createMenu(action, self._SaturneActionIdMap[SaturneHelpUserManualMenu])
        action_id = sgPyQt.actionId(action)
        self._ActionMap[action_id] = action
        self._SaturneActionIdMap[SaturneHelpCodeSaturneAction] = action_id
        self.connect(action, SIGNAL("activated()"), self.slotSaturneHelpCS)

        action = sgPyQt.createAction(SaturneHelpSolutionDomainAction,\
                                      ObjectTR.tr("SATURNE_HELP_SD_ACTION_TEXT"),\
                                      ObjectTR.tr("SATURNE_HELP_SD_ACTION_TIP"),\
                                      ObjectTR.tr("SATURNE_HELP_SD_ACTION_SB"))
        sgPyQt.createMenu(action, self._SaturneActionIdMap[SaturneHelpUserManualMenu])
        action_id = sgPyQt.actionId(action)
        self._ActionMap[action_id] = action
        self._SaturneActionIdMap[SaturneHelpSolutionDomainAction] = action_id
        self.connect(action, SIGNAL("activated()"), self.slotSaturneHelpSD)

        action = sgPyQt.createAction(SaturneHelpCS_KernelAction,\
                                      ObjectTR.tr("SATURNE_HELPCS_KERNEL_ACTION_TEXT"),\
                                      ObjectTR.tr("SATURNE_HELPCS_KERNEL_ACTION_TIP"),\
                                      ObjectTR.tr("SATURNE_HELPCS_KERNEL_ACTION_SB"))
        sgPyQt.createMenu(action, self._SaturneActionIdMap[SaturneHelpUserManualMenu])
        action_id = sgPyQt.actionId(action)
        self._ActionMap[action_id] = action
        self._SaturneActionIdMap[SaturneHelpCS_KernelAction] = action_id
        self.connect(action, SIGNAL("activated()"), self.slotSaturneHelpCS_Kernel)

        action = sgPyQt.createAction(SaturneHelpCS_InfosAction,\
                                      ObjectTR.tr("SATURNE_HELPCS_INFOS_ACTION_TEXT"),\
                                      ObjectTR.tr("SATURNE_HELPCS_INFOS_ACTION_TIP"),\
                                      ObjectTR.tr("SATURNE_HELPCS_INFOS_ACTION_SB"))
        sgPyQt.createMenu(action, self._SaturneActionIdMap[SaturneHelpUserManualMenu])
        action_id = sgPyQt.actionId(action)
        self._ActionMap[action_id] = action
        self._SaturneActionIdMap[SaturneHelpCS_InfosAction] = action_id
        self.connect(action, SIGNAL("activated()"), self.slotSaturneHelpCS_Infos)

        # Neptune Actions
        # Window menu
        action_id = sgPyQt.createMenu(ObjectTR.tr("NEPTUNE_WIN_MENU_TEXT"), 6)
        self._NeptuneActionIdMap[NeptuneWinMenu] = action_id
        action = sgPyQt.createAction(-1,\
                                      ObjectTR.tr("NEPTUNE_WINBROWSER_ACTION_TEXT"),\
                                      ObjectTR.tr("NEPTUNE_WINBROWSER_ACTION_TIP"),\
                                      ObjectTR.tr("NEPTUNE_WINBROWSER_ACTION_SB"))

        action.setChecked(True)  # QAction class
        action.setEnabled(True)

        sgPyQt.createMenu(action, self._NeptuneActionIdMap[NeptuneWinMenu], 100)
        action_id = sgPyQt.actionId(action)
        self._ActionMap[action_id] = action
        self._NeptuneActionIdMap[NeptuneWinBrowserAction] = action_id
        self.connect(action, SIGNAL("toggled(bool)"), self.slotNeptuneWinBrowser)

        action = sgPyQt.createAction(-1,\
                                      ObjectTR.tr("NEPTUNE_WINIDENTITY_ACTION_TEXT"),\
                                      ObjectTR.tr("NEPTUNE_WINIDENTITY_ACTION_TIP"),\
                                      ObjectTR.tr("NEPTUNE_WINIDENTITY_ACTION_SB"))
        action.setChecked(True)
        action.setEnabled(True)
        sgPyQt.createMenu(action, self._NeptuneActionIdMap[NeptuneWinMenu], 100)
        action_id = sgPyQt.actionId(action)
        self._ActionMap[action_id] = action
        self._NeptuneActionIdMap[NeptuneWinIdentityAction] = action_id
        self.connect(action, SIGNAL("toggled(bool)"), self.slotNeptuneWinIdenty)

        action = sgPyQt.createAction(-1,\
                                      ObjectTR.tr("STOP_SOLVER_ACTION_TEXT"),\
                                      ObjectTR.tr("STOP_SOLVER_ACTION_TIP"),\
                                      ObjectTR.tr("STOP_SOLVER_ACTION_SB"))
        action_id = sgPyQt.actionId(action)
        self._ActionMap[action_id] = action
        self._SolverActionIdMap[StopSolverAction] = action_id
        self.connect(action, SIGNAL("activated()"), self.slotStopSolver)

        action = sgPyQt.createAction(-1,\
                                      ObjectTR.tr("SHOW_SOLVER_PROCESS_ACTION_TEXT"),\
                                      ObjectTR.tr("SHOW_SOLVER_PROCESS_ACTION_TIP"),\
                                      ObjectTR.tr("SHOW_SOLVER_PROCESS_ACTION_SB"))
        action_id = sgPyQt.actionId(action)
        self._ActionMap[action_id] = action
        self._SolverActionIdMap[ShowSolverProcessAction] = action_id
        self.connect(action, SIGNAL("activated()"), self.slotShowSolverProcess)

#        action_id = sgPyQt.createMenu(ObjectTR.tr("MESH_OR_GROUP_REPRESENTATION"), -1, -1)
#        self._CommonActionIdMap[SaturneHelpUserManualMenu] = action_id

        action = sgPyQt.createAction(-1,\
                                      ObjectTR.tr("MESH_OR_GROUP_REPRESENTATION_SHADED"),\
                                      ObjectTR.tr("MESH_OR_GROUP_REPRESENTATION_SHADED"),\
                                      ObjectTR.tr("MESH_OR_GROUP_REPRESENTATION_SHADED"))
        action_id = sgPyQt.actionId(action)
        self._ActionMap[action_id] = action
        self._CommonActionIdMap[DisplayTypeSHADED] = action_id
        self.connect(action, SIGNAL("activated()"), self.slotDisplayTypeSHADED)

        action = sgPyQt.createAction(-1,\
                                      ObjectTR.tr("MESH_OR_GROUP_REPRESENTATION_WIREFRAME"),\
                                      ObjectTR.tr("MESH_OR_GROUP_REPRESENTATION_WIREFRAME"),\
                                      ObjectTR.tr("MESH_OR_GROUP_REPRESENTATION_WIREFRAME"))
        action_id = sgPyQt.actionId(action)
        self._ActionMap[action_id] = action
        self._CommonActionIdMap[DisplayTypeWIREFRAME] = action_id
        self.connect(action, SIGNAL("activated()"), self.slotDisplayTypeWIREFRAME)


    def updateActions(self):
        """
        Updates all action according with current selection and study states.
        This function connected to selection change signal.
        """
        if Trace(): print "calling updateActions..."
        component = CFDSTUDYGUI_DataModel._getComponent()
        if component == None:
            #disable all actions except Study Location
            for i in self._CommonActionIdMap:
                if not i == InfoCFDSTUDYAction:
                    if i == SetStudyAction:
                        self.commonAction(i).setEnabled(True)
                    else:
                        self.commonAction(i).setEnabled(False)
        else:
            #enable all actions
            for i in self._CommonActionIdMap:
                if not i == InfoCFDSTUDYAction:
                    self.commonAction(i).setEnabled(True)

        # selection handler
        sobj = self._singleSelectedObject()
        isStudy = CFDSTUDYGUI_DataModel.checkType(sobj, CFDSTUDYGUI_DataModel.dict_object["Study"])
        self.commonAction(AddCaseAction).setEnabled(isStudy)
        self.commonAction(RunCaseAction).setEnabled(isStudy)
        aStudy = CFDSTUDYGUI_DataModel.GetStudyByObj(sobj)
        aCase = CFDSTUDYGUI_DataModel.GetCase(sobj)

        if aStudy != None and aCase != None:
            self.commonAction(LaunchGUIAction).setEnabled(CFDSTUDYGUI_DataModel.checkCaseLaunchGUI(aCase))
            self.commonAction(OpenXMLCFDGUIAction).setEnabled(CFDSTUDYGUI_DataModel.checkCaseLaunchGUI(aCase))
        else:
            self.commonAction(LaunchGUIAction).setEnabled(False)

        #enable / disable solver actions
        #isActivatedView = True # temp solution
        isActivatedView = self._SolverGUI.isActive() # Main GUI Window is active

        for a in self._SolverActionIdMap:
            if a != SolverFileMenu and a != SolverToolsMenu \
                and a != SolverHelpMenu and a != StopSolverAction \
                and a != ShowSolverProcessAction:
                self.solverAction(a).setEnabled(isActivatedView)

        if CFD_Code() == CFD_Saturne:
            for a in self._SaturneActionIdMap:
                if a != SaturneHelpUserManualMenu:
                    self.solverAction(a).setEnabled(isActivatedView)
        elif CFD_Code() == CFD_Neptune:
            for a in self._NeptuneActionIdMap:
                if a != NeptuneWinMenu:
                    self.solverAction(a).setEnabled(isActivatedView)


    def customPopup(self, id, popup):
        """
        Callback for fill popup menu according current selection state.

        @type id: C{int}
        @param id: type of the branch tree slected in the Object Brower.
        @type popup: C{QPopupMenu}
        @param popup: popup menu from the Object Browser.
        """
        if id == CFDSTUDYGUI_DataModel.dict_object["Study"]:
            popup.addAction(self.commonAction(AddCaseAction))
            popup.addAction(self.commonAction(UpdateObjBrowserAction))
            popup.addAction(self.commonAction(RunCaseAction))
        elif id == CFDSTUDYGUI_DataModel.dict_object["Case"]:
            popup.addAction(self.commonAction(LaunchGUIAction))
            popup.addAction(self.commonAction(RemoveAction))
        elif id == CFDSTUDYGUI_DataModel.dict_object["PRETLink"]:
            popup.addAction(self.commonAction(DefineLinkAction))
            popup.addAction(self.commonAction(RemoveAction))
            popup.addAction(self.commonAction(UpdateObjBrowserAction))
            popup.addAction(self.commonAction(LaunchGUIAction))
        elif id == CFDSTUDYGUI_DataModel.dict_object["SUITELink"]:
            popup.addAction(self.commonAction(DefineLinkAction))
            popup.addAction(self.commonAction(RemoveAction))
        elif id == CFDSTUDYGUI_DataModel.dict_object["DATAFile"]:
            popup.addAction(self.commonAction(EditAction))
            popup.addAction(self.commonAction(MoveToDRAFTAction))
            popup.addAction(self.commonAction(CopyCaseFileAction))
        elif id == CFDSTUDYGUI_DataModel.dict_object["DATADRAFTFile"]:
            popup.addAction(self.commonAction(EditAction))
            popup.addAction(self.commonAction(RemoveAction))
            popup.addAction(self.commonAction(CopyInDATAAction))
        elif id == CFDSTUDYGUI_DataModel.dict_object["THCHFile"]:
            popup.addAction(self.commonAction(ViewAction))
            popup.addAction(self.commonAction(CopyInDATAAction))
            popup.addAction(self.commonAction(CopyCaseFileAction))
        elif id == CFDSTUDYGUI_DataModel.dict_object["DATALaunch"]:
            # popup associated with SaturneGUI  script file under DATA directory
            popup.addAction(self.commonAction(LaunchGUIAction))
        elif id == CFDSTUDYGUI_DataModel.dict_object["DATAfileXML"]:
            popup.addAction(self.commonAction(OpenXMLCFDGUIAction))
            popup.addAction(self.commonAction(CopyCaseFileAction))
        elif id == CFDSTUDYGUI_DataModel.dict_object["SRCFolder"]:
            popup.addAction(self.commonAction(CheckCompilationAction))
        elif id == CFDSTUDYGUI_DataModel.dict_object["SRCFile"]:
            popup.addAction(self.commonAction(EditAction))
            popup.addAction(self.commonAction(MoveToDRAFTAction))
            popup.addAction(self.commonAction(CopyCaseFileAction))
        elif id == CFDSTUDYGUI_DataModel.dict_object["SRCDRAFTFile"]:
            popup.addAction(self.commonAction(EditAction))
            popup.addAction(self.commonAction(RemoveAction))
            popup.addAction(self.commonAction(CopyInSRCAction))
        elif id == CFDSTUDYGUI_DataModel.dict_object["LOGSRCFile"]:
            popup.addAction(self.commonAction(ViewAction))
        elif id == CFDSTUDYGUI_DataModel.dict_object["USRSRCFile"]:
            popup.addAction(self.commonAction(ViewAction))
            popup.addAction(self.commonAction(CopyInSRCAction))
        elif id == CFDSTUDYGUI_DataModel.dict_object["RESUFile"]:
            popup.addAction(self.commonAction(ViewAction))
        elif id == CFDSTUDYGUI_DataModel.dict_object["VirtFolder"]:
            popup.addAction(self.commonAction(RemoveAction))
        elif id == CFDSTUDYGUI_DataModel.dict_object["RESSRCFile"]:
            popup.addAction(self.commonAction(ViewAction))
        elif id == CFDSTUDYGUI_DataModel.dict_object["HISTFile"]:
            popup.addAction(self.commonAction(ViewAction))
            popup.addAction(self.commonAction(ExportInPostProAction))
        elif id == CFDSTUDYGUI_DataModel.dict_object["PRETFolder"] or \
             id == CFDSTUDYGUI_DataModel.dict_object["SUITEFolder"]:
            popup.addAction(self.commonAction(CreateLinkAction))
            popup.addAction(self.commonAction(CopyInDATAAction))
        elif id == CFDSTUDYGUI_DataModel.dict_object["RESMEDFile"]:
            popup.addAction(self.commonAction(ExportInPostProAction))
        elif id == CFDSTUDYGUI_DataModel.dict_object["SCRPTLanceFile"]:
            popup.addAction(self.commonAction(ViewAction))
            popup.addAction(self.commonAction(RunScriptAction))
        elif id == CFDSTUDYGUI_DataModel.dict_object["SCRPTScriptFile"]:
            popup.addAction(self.commonAction(EditAction))
            popup.addAction(self.commonAction(RunScriptAction))
        elif id == CFDSTUDYGUI_DataModel.dict_object["SCRPTFile"]:
            popup.addAction(self.commonAction(ViewAction))
        elif id == CFDSTUDYGUI_DataModel.dict_object["SCRPTStdLog"]:
            popup.addAction(self.solverAction(StopSolverAction))
        elif id == CFDSTUDYGUI_DataModel.dict_object["SCRPTExtLog"]:
            popup.addAction(self.solverAction(StopSolverAction))
            popup.addAction(self.solverAction(ShowSolverProcessAction))
        elif id == CFDSTUDYGUI_DataModel.dict_object["FICHEFile"]:
            popup.addAction(self.commonAction(ViewAction))
        elif id == CFDSTUDYGUI_DataModel.dict_object["DESFile"] \
             or id == CFDSTUDYGUI_DataModel.dict_object["CGNSFile"] \
             or id == CFDSTUDYGUI_DataModel.dict_object["GeomFile"] \
             or id == CFDSTUDYGUI_DataModel.dict_object["CaseFile"] \
             or id == CFDSTUDYGUI_DataModel.dict_object["NeuFile"] \
             or id == CFDSTUDYGUI_DataModel.dict_object["MSHFile"] \
             or id == CFDSTUDYGUI_DataModel.dict_object["HexFile"] \
             or id == CFDSTUDYGUI_DataModel.dict_object["UnvFile"]:
            popup.addAction(self.commonAction(ECSConvertAction))
        elif id == CFDSTUDYGUI_DataModel.dict_object["MEDFile"]:
            popup.addAction(self.commonAction(ExportInSMESHAction))
        elif id == CFDSTUDYGUI_DataModel.dict_object["MESHFile"]:
            popup.addAction(self.commonAction(ViewAction))
        elif id == CFDSTUDYGUI_DataModel.dict_object["DATFile"]:
            popup.addAction(self.commonAction(EditAction))
            #popup.addAction(self.commonAction(ECSConvertAction))
        elif id == CFDSTUDYGUI_DataModel.dict_object["POSTFile"]:
            popup.addAction(self.commonAction(ViewAction))
        elif id == CFDSTUDYGUI_DataModel.dict_object["PublishedMeshIntoObjectBrowser"]:
            popup.addAction(self.commonAction(DisplayMESHAction))
            #popup.addAction(self.commonAction(DisplayOnlyMESHAction))
            #popup.addAction(self.commonAction(HideMESHAction))
        elif id == CFDSTUDYGUI_DataModel.dict_object["PublishedMeshGroupIntoObjectBrowser"]:
            popup.addAction(self.commonAction(DisplayGroupMESHAction))
            popup.addAction(self.commonAction(DisplayOnlyGroupMESHAction))
            popup.addAction(self.commonAction(HideGroupMESHAction))

        elif id == "VTKViewer":
            popup.addAction(self.commonAction(DisplayTypeSHADED))
            popup.addAction(self.commonAction(DisplayTypeWIREFRAME))


    def slotStudyLocation(self):
        """
        Loads the CFD study location. If the name of the CFD study
        does not exists, the corresponding folder is created.
        """
        if Trace(): print 'CFDSTUDYGUI_ActionsHandler.slotStudyLocation'

        dialog = self.DialogCollector.SetTreeLocationDialog
        dialog.exec_()
        if not self.DialogCollector.SetTreeLocationDialog.result() == QDialog.Accepted:
            return

        cursor = QCursor(Qt.BusyCursor)
        QApplication.setOverrideCursor(cursor)

        self._CommandMgr.runFunctionDlg(CFDSTUDYGUI_DataModel._SetStudyLocation,
                                        self.tr("STMSG_SET_STUDY_LOCATION"),\
                                        False,\
                                        theStudyPath = dialog.StudyPath, \
                                        theCaseNames = dialog.CaseNames)
        studyId = sgPyQt.getStudyId()
        sgPyQt.updateObjBrowser(studyId, 1)
        self.updateActions()
        # self.enableHelpSaturne()
        QApplication.restoreOverrideCursor()


    def enableHelpSaturne(self):
        """
        """
        #FIXME: a implementer avec les tr
        for a in self._SaturneActionIdMap:
            if self.solverAction(a).text() in ["About",
                                               "License",
                                               "Code_Saturne",
                                               "Solution domain",
                                               "Code_Saturne kernel",
                                               "Code_Saturne infos"]:
                self.solverAction(a).setEnabled(True)


    def slotAddCase(self):
        """
        Builds new CFD cases.
        """
        dialog = self.DialogCollector.SetTreeLocationDialog
        dialog.setCaseMode(True)

        studyObj = self._singleSelectedObject()
        if studyObj == None:
            return

        dialog.StudyPath = CFDSTUDYGUI_DataModel._GetPath(studyObj)
        dialog.exec_()
        if  not self.DialogCollector.SetTreeLocationDialog.result() == QDialog.Accepted:
            #cancel of new case creation
            dialog.setCaseMode(False)
            return

        dialog.setCaseMode(False)

        cursor = QCursor(Qt.BusyCursor)
        QApplication.setOverrideCursor(cursor)

        self._CommandMgr.runFunctionDlg(CFDSTUDYGUI_DataModel._SetStudyLocation,
                                        self.tr("STMSG_ADD_CASE"),\
                                        False,\
                                        theStudyPath = dialog.StudyPath, \
                                        theCaseNames = dialog.CaseNames)

        QApplication.restoreOverrideCursor()


    def slotInfo(self):
        """
        Shows the QDialog with the info from CFDSTUDY:
            - CFD code selected
            - environnement variables defined
        """
        if Trace(): print "CFDSTUDYGUI_ActionsHandler.slotInfo"
        dialog = self.DialogCollector.InfoDialog
        dialog.show()
        self.updateActions()


    def slotUpdateObjectBrowser(self):
        """
        Re-reads the unix folders and updates the complete representation
        of the CFD studies in the Object Browser.
        """
        self.updateObjBrowser()


    def updateObjBrowser(self, Object=None):
        """
        Updates CFD study sub-tree from the argument object.

        @type theObject: C{SObject}
        @param theObject: branch of a tree of data to update.
        """
        if Trace(): print "CFDSTUDYGUI_ActionsHandler.updateObjBrowser"

        self._CommandMgr.runFunctionDlg(CFDSTUDYGUI_DataModel.UpdateSubTree,\
                                         self.tr("STMSG_UPDATE_STUDY"),\
                                         False,\
                                         theObject = Object)

    def slotViewAction(self):
        """
        Edits in the read only mode the file selected in the Object Browser.
        Warning, the editor is always emacs!
        """

        viewer = self.tr("CFDSTUDY_PREF_READER")
        if not viewer.isEmpty():
            viewerName = str(viewer.toLatin1())
            sobj = self._singleSelectedObject()
            if sobj is not None:
                path = CFDSTUDYGUI_DataModel._GetPath(sobj)
                if Trace(): print "PATH:", path
                if re.match(".*emacs$", viewerName):
                    if Trace(): print "Read Only"
                    os.spawnlp(os.P_NOWAIT,viewerName , viewerName, path, "-f", "toggle-read-only")
                else:
                    os.spawnlp(os.P_NOWAIT,viewerName ,viewerName , path)


    def slotEditAction(self):
        """
        Edits in the user's editor the file selected in the Object Browser.
        """
        viewer = self.tr("CFDSTUDY_PREF_EDITOR")

        if not viewer.isEmpty():
            viewerName = str(viewer.toLatin1())
            sobj = self._singleSelectedObject()
            if not sobj == None:
                path = CFDSTUDYGUI_DataModel._GetPath(sobj)
                os.spawnlp(os.P_NOWAIT,viewerName ,viewerName , path)


    def slotRemoveAction(self):
        """
        Deletes file or folder from the Object Browser, and from the unix system files.
        Delete dock windows attached to a CFD Study if this study is deleted from the Object Browser.
        """
        sobj = self._singleSelectedObject()
        if not sobj == None:
            mess = ObjectTR.tr("REMOVE_ACTION_CONFIRM_MESS").arg(sobj.GetName())
            if QMessageBox.warning(None, "Warning", mess, QMessageBox.Yes, QMessageBox.No) == QMessageBox.No:
                return

            watchCursor = QCursor(Qt.WaitCursor)
            QApplication.setOverrideCursor(watchCursor)

            path = CFDSTUDYGUI_DataModel._GetPath(sobj)

            caseName = sobj.GetName()
            father = sobj.GetFather()
            fatherpath = CFDSTUDYGUI_DataModel._GetPath(father)
            fathername = father.GetName()
            if Trace(): print 'CFDSTUDYGUI_ActionHandler.slotRemoveAction', father, fatherpath, fathername
            if fathername == 'RESU':
                name = '*' + sobj.GetName()
                path = os.path.join(fatherpath , name)
            rmpath = 'rm -fr ' + path
            if Trace(): print 'slotRemoveAction rm -fr',path
            os.system(rmpath)
            if os.listdir(fatherpath) == []:
                rmfatherpath = 'rm -fr ' + fatherpath
                if Trace(): print 'slotRemoveAction rm -fr',fatherpath

            CFDSTUDYGUI_DataModel._RebuildTreeRecursively(sobj)
            CFDSTUDYGUI_SolverGUI.updateObjectBrowser()
            CFDSTUDYGUI_SolverGUI.removeDockWindow(caseName)
            QApplication.restoreOverrideCursor()


    def slotCopyInDATA(self):
        """
        """
        sobj = self._singleSelectedObject()
        if sobj != None:
            studyId = sgPyQt.getStudyId()
            path = CFDSTUDYGUI_DataModel._GetPath(sobj)

            study = CFDSTUDYGUI_DataModel._getStudy()
            builder = study.NewBuilder()

            attr = builder.FindOrCreateAttribute(sobj, "AttributeLocalID")
            if attr.Value() == CFDSTUDYGUI_DataModel.dict_object["PRETFolder"] or \
               attr.Value() == CFDSTUDYGUI_DataModel.dict_object["SUITEFolder"]:
                case = CFDSTUDYGUI_DataModel.GetCase(sobj)
                if not case == None:
                    iter  = study.NewChildIterator(case)
                    while iter.More():
                        if Trace(): print "Currrent Name: ", iter.Value().GetName()
                        if iter.Value().GetName() == "DATA":
                            newpath = os.path.join(CFDSTUDYGUI_DataModel._GetPath(iter.Value()), sobj.GetName())
                            #remove if exists

                            if os.path.exists(newpath):
                                mess = ObjectTR.tr("OVERWRITE_CONFIRM_MESS").arg(sobj.GetName())
                                if QMessageBox.warning(None, "Warning", mess, QMessageBox.Yes, QMessageBox.No) == QMessageBox.No:
                                    return
                                os.spawnlp(os.P_WAIT, 'rm', 'rm', '-fr', newpath)

                            os.spawnlp(os.P_WAIT, 'cp', 'cp', '-fR', path, newpath)
                            if Trace(): print "coping is successfull!!!"

                            #ubdate Object Browser
                            CFDSTUDYGUI_DataModel._RebuildTreeRecursively(case)
                            break
                        iter.Next()
                    sgPyQt.updateObjBrowser(studyId, 1)
            else:
                parent = sobj.GetFather()
                if not parent == None:
                    parent = parent.GetFather()
                    if not parent == None and parent.GetName() == "DATA":
                        parentPath = CFDSTUDYGUI_DataModel._GetPath(parent)
                        newpath = os.path.join(parentPath , sobj.GetName())

                        if os.path.exists(newpath):
                            mess = ObjectTR.tr("OVERWRITE_CONFIRM_MESS")
                            if QMessageBox.warning(None, "Warning", mess, QMessageBox.Yes, QMessageBox.No) == QMessageBox.No:
                                return
                        os.spawnlp(os.P_WAIT, 'cp', 'cp', '-f', path, parentPath)

                        #ubdate Object Browser
                        CFDSTUDYGUI_DataModel._RebuildTreeRecursively(parent)
                        sgPyQt.updateObjBrowser(studyId, 1)
                        if Trace(): print "coping is successfull!!!"
                    else:
                        if Trace(): print "No DATA parent"
                else:
                    if Trace(): print "No parent"


    def slotCopyInSRC(self):
        """
        """
        sobj = self._singleSelectedObject()
        if sobj != None:
            path = CFDSTUDYGUI_DataModel._GetPath(sobj)
            parent = sobj.GetFather()
            if not parent == None:
                if not parent == None and (parent.GetName() != "DRAFT" and parent.GetName() != "REFERENCE"):
                    parent = parent.GetFather()
                if not parent == None and (parent.GetName() == "REFERENCE" or parent.GetName() == "DRAFT"):
                    parent = parent.GetFather()
                    if not parent == None and parent.GetName() == "SRC":
                        parentPath = CFDSTUDYGUI_DataModel._GetPath(parent)
                        destPath = parentPath + '/' + sobj.GetName()
                        if os.path.exists(destPath):
                            mess = ObjectTR.tr("OVERWRITE_CONFIRM_MESS")
                            if QMessageBox.warning(None, "Warning", mess, QMessageBox.Yes, QMessageBox.No) == QMessageBox.No:
                                return
                        os.spawnlp(os.P_WAIT, 'cp', 'cp', '-f', path, parentPath)

                        #update Object Browser
                        CFDSTUDYGUI_DataModel._RebuildTreeRecursively(parent)
                        studyId = sgPyQt.getStudyId()
                        sgPyQt.updateObjBrowser(studyId, 1)
                    else:
                        if Trace(): print "No SRC parent"
                else:
                    if Trace(): print "No USERS or DRAFT parent"
            else:
                if Trace(): print "No parent"


    def slotMoveToDRAFT(self):
        """
        """
        sobj = self._singleSelectedObject()
        if sobj != None:
            path = CFDSTUDYGUI_DataModel._GetPath(sobj)
            parent = sobj.GetFather()
            if not parent == None:
                parentPath = CFDSTUDYGUI_DataModel._GetPath(parent) + '/DRAFT'
                destPath = parentPath + '/' + sobj.GetName()
                if os.path.exists(destPath):
                    mess = ObjectTR.tr("OVERWRITE_CONFIRM_MESS")
                    if QMessageBox.warning(None, "Warning", mess, QMessageBox.Yes, QMessageBox.No) == QMessageBox.No:
                        return

                if os.path.exists(parentPath) == False:
                    os.mkdir(parentPath)

                os.spawnlp(os.P_WAIT, 'mv', 'mv', '-f', path, parentPath)
                CFDSTUDYGUI_DataModel._RebuildTreeRecursively(parent)
                studyId = sgPyQt.getStudyId()
                sgPyQt.updateObjBrowser(studyId, 1)
            else:
                if Trace(): print "No parent"


    def _singleSelectedObject(self):
        """
        """
        if Trace(): "selected count ", sg.SelectedCount()
        study = CFDSTUDYGUI_DataModel._getStudy()
        if sg.SelectedCount() == 1:
            entry = sg.getSelected(0)
            if entry != '':
                if Trace(): print "Entry:", entry
                return study.FindObjectID(entry)

        return None


    def _multipleSelectedObject(self):
        """
        """
        if Trace(): "selected count ", sg.SelectedCount()
        study = CFDSTUDYGUI_DataModel._getStudy()

        i = 0
        liste_SObj = []
        while i < sg.SelectedCount():
            entry = sg.getSelected(i)
            if entry != '':
                if Trace():
                    print "Entry:", entry
                liste_SObj.append(study.FindObjectID(entry))
            i = i+1
        return liste_SObj


    def slotExportInPostPro(self):
        """
        """
        waitCursor = QCursor(Qt.WaitCursor)
        QApplication.setOverrideCursor(waitCursor)

        sobj = self._singleSelectedObject()
        if sobj != None:
            path = CFDSTUDYGUI_DataModel._GetPath(sobj)
            if re.match(".*\.med$", sobj.GetName()):
                #export Med file
                self.myVisu.ImportFile(path)
            elif re.match(".*\.dat$", sobj.GetName()):
                self.myVisu.ImportTables(path)
            studyId = sgPyQt.getStudyId()
            sgPyQt.updateObjBrowser(studyId, 1)
        QApplication.restoreOverrideCursor()


    def slotExportInSMESH(self):
        """
        """
        waitCursor = QCursor(Qt.WaitCursor)
        QApplication.setOverrideCursor(waitCursor)
        from smesh import FACE, VOLUME

        sobj = self._singleSelectedObject()
        if sobj != None:
            path = CFDSTUDYGUI_DataModel._GetPath(sobj)
            smesh_component = salome.lcc.FindOrLoadComponent("FactoryServer", "SMESH")
            smesh_component.SetCurrentStudy(salome.myStudy)

            if smesh_component and re.match(".*\.med$", sobj.GetName()):
                #export Med file
                aMeshes, aStatus = smesh_component.CreateMeshesFromMED(path)
                if not aStatus:
                    QApplication.restoreOverrideCursor()
                    mess = ObjectTR("EXPORT_IN_SMESH_ACTION_WARNING")
                    QMessageBox.warning(None, "Warning", mess, QMessageBox.Ok, 0)
                    return

#                SO_father = salome.myStudy.FindComponent("CFDSTUDY")
#                if not SO_father == None:
#                    for mesh in aMeshes:
#                        idmesh = CFDSTUDYGUI_DataModel.dict_object["PublishedMeshIntoObjectBrowser"]
#                        SO_mesh = CFDSTUDYGUI_DataModel.publishInStudySalome(sobj, mesh.GetName(), idmesh)
#
#                        for group in mesh.GetGroups():
#                            if group.GetType() == FACE:
#                                id = CFDSTUDYGUI_DataModel.dict_object["PublishedMeshGroupFacesIntoObjectBrowser"]
#                                SO_gfaces = CFDSTUDYGUI_DataModel.publishInStudySalome(SO_mesh, "Groups of Faces", id)
#                                break
#                        for group in mesh.GetGroups():
#                            if group.GetType() == VOLUME:
#                                id = CFDSTUDYGUI_DataModel.dict_object["PublishedMeshGroupCellsIntoObjectBrowser"]
#                                SO_gcells = CFDSTUDYGUI_DataModel.publishInStudySalome(SO_mesh, "Groups of Cells", id)
#                                break
#
#                        for group in mesh.GetGroups():
#                            if Trace(): print "slotExportInSMESH: group.GetName() = ",group.GetName()
#                            if group.GetType() == FACE:
#                                idgroup = CFDSTUDYGUI_DataModel.dict_object["PublishedMeshGroupIntoObjectBrowser"]
#                                CFDSTUDYGUI_DataModel.publishInStudySalome(SO_gfaces, group.GetName(), idgroup)
#                            if group.GetType() == VOLUME:
#                                idgroup = CFDSTUDYGUI_DataModel.dict_object["PublishedMeshGroupIntoObjectBrowser"]
#                                CFDSTUDYGUI_DataModel.publishInStudySalome(SO_gcells, group.GetName(), idgroup)
#
            studyId = sgPyQt.getStudyId()
            sgPyQt.updateObjBrowser(studyId, 1)

        QApplication.restoreOverrideCursor()


    def slotDisplayMESH(self):
        """
        """
        waitCursor = QCursor(Qt.WaitCursor)
        QApplication.setOverrideCursor(waitCursor)

        self.myView = self.myViewManager.GetCurrentView()

        sobj = self._singleSelectedObject()
        if sobj != None:
            father = sobj.GetFather()
            medFilePath = CFDSTUDYGUI_DataModel._GetPath(father)
            if re.match(".*\.med$", father.GetName()):
                myResult = self.myVisu.ImportFile(medFilePath)
                if myResult is None:
                    print "Error"
                    return
                mesh = self.myVisu.MeshOnEntity(myResult, str(sobj.GetName()), VISU.CELL)
                aColor = SALOMEDS.Color(0,1,1)
                mesh.SetCellColor(aColor)
                mesh.SetPresentationType(VISU.SURFACEFRAME)
                self.myView.Display(mesh)

                self.myView.Update()
                self.myView.FitAll()
                studyId = sgPyQt.getStudyId()
                sgPyQt.updateObjBrowser(studyId, 1)

        QApplication.restoreOverrideCursor()


    def slotDisplayMESH2(self):
        """
        """
        waitCursor = QCursor(Qt.WaitCursor)
        QApplication.setOverrideCursor(waitCursor)

        import smesh
        smeshgui = salome.ImportComponentGUI("SMESH")
        smeshgui.Init(salome.myStudyId)

        study = CFDSTUDYGUI_DataModel._getStudy()
        builder = study.NewBuilder()

        sobj = self._singleSelectedObject()
        if not sobj == None:
            attr = builder.FindOrCreateAttribute(sobj, "AttributeLocalID")
            if attr.Value() == CFDSTUDYGUI_DataModel.dict_object["PublishedMeshIntoObjectBrowser"]:
                log.debug("slotDisplayMESH -> entry = %s" % sobj.GetID())
                smeshgui.CreateAndDisplayActor(sobj.GetID())
                #smeshgui.CreateAndDisplayActor('0:1:3:3')
                #salome.sg.Display(sobj.GetID())
                smesh.UpdateView()
                salome.sg.UpdateView()
                salome.sg.FitAll()

        QApplication.restoreOverrideCursor()


    def slotDisplayMESHGroups(self):
        self.__displayMeshGroups("display", VISU.SURFACEFRAME)


    def slotDisplayOnlyMESHGroups(self):
        self.__displayMeshGroups("only", VISU.SURFACEFRAME)


    def slotHideMESHGroups(self):
        self.__displayMeshGroups("hide", VISU.SURFACEFRAME)


    def __displayMeshGroups(self, display, type=None):
        """
        """
        waitCursor = QCursor(Qt.WaitCursor)
        QApplication.setOverrideCursor(waitCursor)

        self.myView = self.myViewManager.GetCurrentView()

        for group in self._multipleSelectedObject():
            group_type = group.GetFather()
            fatherMesh = group_type.GetFather()
            fatherFile = fatherMesh.GetFather()
            medFilePath = CFDSTUDYGUI_DataModel._GetPath(fatherFile)
            if re.match(".*\.med$", fatherFile.GetName()):
                myResult = self.myVisu.ImportFile(medFilePath)
                if myResult is None:
                    print "slotDisplayMESHGroups: Error: myResult is None"
                    return
                mesh_group = self.myVisu.GroupMesh(myResult,
                                                   str(fatherMesh.GetName()),
                                                   str(group.GetName()))
                mesh_group.SetPresentationType(type)
                if display == "display":
                    pass
                    #self.myView.Display(mesh_group)
                elif display == "only":
                    self.myView.EraseAll()
                    self.myView.DisplayOnly(mesh_group)
                elif display == "hide":
                    self.myView.Erase(mesh_group)

            self.myView.Update()
            self.myView.FitAll()
        studyId = sgPyQt.getStudyId()
        sgPyQt.updateObjBrowser(studyId, 1)

        QApplication.restoreOverrideCursor()

# VISU.POINT, VISU.WIREFRAME, VISU.SHADED, VISU.INSIDEFRAME,
# VISU.SURFACEFRAME, VISU.FEATURE_EDGES, VISU.SHRINK

    def slotDisplayTypeWIREFRAME(self):
        self.__displayMeshGroups("display", VISU.WIREFRAME)


    def slotDisplayTypeSHADED(self):
        self.__displayMeshGroups("display", VISU.SHADED)


    def slotCreateLink(self):
        """
        """
        sobj = self._singleSelectedObject()
        if sobj != None:
            path = CFDSTUDYGUI_DataModel._GetPath(sobj)
            study = CFDSTUDYGUI_DataModel._getStudy()
            case = CFDSTUDYGUI_DataModel.GetCase(sobj)
            if not case == None:
                iter  = study.NewChildIterator(case)
                while iter.More():
                    if Trace(): print "Currrent Name: ", iter.Value().GetName()
                    if iter.Value().GetName() == "DATA":
                        dataPath = CFDSTUDYGUI_DataModel._GetPath(iter.Value())
                        linkPath = os.path.join(dataPath , sobj.GetName())
                        if Trace(): print "Link Path: ", linkPath
                        if os.path.exists(linkPath):
                            if Trace(): print 'slotCreateLink',linkPath
                            os.spawnlp(os.P_WAIT, 'rm', 'rm', '-fr', linkPath)

                        os.symlink(path, linkPath)
                        #update Object Browser
                        CFDSTUDYGUI_DataModel._RebuildTreeRecursively(case)
                        sgPyQt.updateObjBrowser(1)
                        #CFDSTUDYGUI_SolverGUI.updateObjectBrowser()
                        return
                iter.Next()
            studyId = sgPyQt.getStudyId()
            sgPyQt.updateObjBrowser(studyId, 1)


    def slotDefineLink(self):
        sobj = self._singleSelectedObject()
        if not sobj == None:
            #detection kind of link
            link_kind = -1
            if sobj.GetName() == "PRE_TRAITEMENT":
                link_kind = CFDSTUDYGUI_DataModel.dict_object["PRETLink"]
            elif sobj.GetName() == "SUITE":
                link_kind = CFDSTUDYGUI_DataModel.dict_object["SUITELink"]
            else:
                return

            path = CFDSTUDYGUI_DataModel._GetPath(sobj)
            study = CFDSTUDYGUI_DataModel._getStudy()
            builder = study.NewBuilder()
            case = CFDSTUDYGUI_DataModel.GetCase(sobj)
            if not case == None:
                iter  = study.NewChildIterator(case)
                while iter.More():
                    if Trace():  print "Currrent Name: ", iter.Value().GetName()
                    if iter.Value().GetName() == "RESU":
                        resupath = CFDSTUDYGUI_DataModel._GetPath(iter.Value())
                        ld = os.listdir(resupath)

                        iter1  = study.NewChildIterator(iter.Value())
                        NameList = []
                        while iter1.More():
                            attr = builder.FindOrCreateAttribute(iter1.Value(), "AttributeLocalID")
                            if attr.Value() == CFDSTUDYGUI_DataModel.dict_object["VirtFolder"] and \
                               ld.count(sobj.GetName() + "." + iter1.Value().GetName()) == 1:
                                NameList.append(iter1.Value().GetName())
                            iter1.Next()

                        self.DialogCollector.DefineLinkDialog.fillDialog(NameList)
                        self.DialogCollector.DefineLinkDialog.exec_()
                        status = self.DialogCollector.DefineLinkDialog.result()
                        name = self.DialogCollector.DefineLinkDialog.currentName()

                        if status == QDialog.Accepted and not name == "":
                            #remove old link
                            os.unlink(path)

                            os.symlink(os.path.join(resupath, sobj.GetName()) + "." + str(name.toLatin1()) , path)
                            #update Object Browser
                            CFDSTUDYGUI_DataModel._RebuildTreeRecursively(case)
                            return
                        else:
                            break
                    iter.Next()
                studyId = sgPyQt.getStudyId()
                sgPyQt.updateObjBrowser(studyId, 1)


    def slotOpenCFD_GUI(self) :
        """
        Open into Salome the CFD GUI from an XML file whose name is sobj.GetName()
        """
        sobj = self._singleSelectedObject()
        if Trace(): print "CFDSTUDYGUI_ActionsHandler : slotOpenCFD_GUI : sobj: ", sobj
        if sobj != None :
            if CFDSTUDYGUI_DataModel.checkType(sobj, CFDSTUDYGUI_DataModel.dict_object["DATAfileXML"]) :
                aXmlFileName = sobj.GetName()
                aCase = CFDSTUDYGUI_DataModel.GetCase(sobj)
                aStudy = CFDSTUDYGUI_DataModel.GetStudyByObj(sobj)
                if aCase != None :
                    aCaseName = aCase.GetName()

                else :
                    mess = "Error : "+ aXmlFileName + " file has no CFD Case into the Salome Object browser ; can't open CFD GUI with this file"
                    QMessageBox.warning(None, "Warning", mess, QMessageBox.Ok, 0)
                    return
                if  aStudy != None :

                    aStudyName = aStudy.GetName()
                else :
                    mess = "Error : "+ aXmlFileName + " file has no CFD Study into the Salome Object browser ; can't open CFD GUI with this file"
                    QMessageBox.warning(None, "Warning", mess, QMessageBox.Ok, 0)
                    return
                if CFDSTUDYGUI_SolverGUI.findDockWindow(aXmlFileName, aCaseName,aStudyName) :
                    mess = aStudyName + " " + aCaseName + " : " + aXmlFileName + " is already opened"
                    QMessageBox.information(None, "Information", mess, QMessageBox.Ok, QMessageBox.NoButton)
                    return
                # xml case file not already opened

                aCmd = []
                aCmd.append('-f')
                if CFD_Code() == CFD_Saturne:
                    aCmd.append(aXmlFileName)
                elif CFD_Code() == CFD_Neptune:
                    import os.path
                    aCmd.append(os.path.join("DATA", aXmlFileName))

                wm = self._SolverGUI.ExecGUI(self.dskAgent().workspace(), aXmlFileName, aCase, aCmd)
                self.updateActions()


    def slotRunGUI(self, study=None, case=None):
        """
        Build the command line for the GUI of Code_Saturne/NEPTUNE_CFD.

        @type study: C{String}
        @param study: name of the study.
        @type case: C{String}
        @param case: name of the case.
        """
        #get current selection
        sobj = self._singleSelectedObject()
        if Trace(): print "slotRunGUI : sobj: ", sobj

        if study:
            aStudy = study
        else:
            aStudy = CFDSTUDYGUI_DataModel.GetStudyByObj(sobj)

        if aStudy == None:
            if Trace(): print "Can't find CFD Study object"
            return
        objname = sobj.GetName()

        # get current case
        if case:
            aCase = case
        else:
            aCase = CFDSTUDYGUI_DataModel.GetCase(sobj)

        if Trace(): print "CFDSTUDYGUI_ActionsHandler::::::::::::::::::: slotRunGUI Case: ", aCase
        self.DialogCollector.GUIActivationDialog.setCurrentCase(aCase)
        self.DialogCollector.GUIActivationDialog.setCurrentStudy(aStudy)
        if CFDSTUDYGUI_DataModel.checkType(sobj, CFDSTUDYGUI_DataModel.dict_object["DATAfileXML"]) :
            xmlFileName = objname
        else :
            xmlFileName = ""
        self.DialogCollector.GUIActivationDialog.fillData(xmlFileName)

        if aCase != None:
            if CFDSTUDYGUI_DataModel.checkType(sobj, CFDSTUDYGUI_DataModel.dict_object["DATAfileXML"]) or \
               CFDSTUDYGUI_DataModel.checkType(sobj, CFDSTUDYGUI_DataModel.dict_object["RESXMLFile"]):
                #checks that such tab already opened

                # launch GUI from an xml file from the Object browser

                if CFDSTUDYGUI_SolverGUI.findDockWindow(sobj.GetName(), aCase.GetName(),aStudy.GetName()) :
                    mess = aStudy.GetName() + " " + aCase.GetName() + " : " + sobj.GetName() + " is already launched"
                    QMessageBox.information(None, "Information", mess, QMessageBox.Ok, QMessageBox.NoButton)
                    return
                # xml case file not already opened
                self.DialogCollector.GUIActivationDialog.setCurrentXMLfile(sobj.GetName())
                self.updateActions()
        self.DialogCollector.GUIActivationDialog.exec_()
        if not self.DialogCollector.GUIActivationDialog.result() == QDialog.Accepted:
            return

        aCaseName = self.DialogCollector.GUIActivationDialog.currentCaseName()
        aXmlFile = None
        aCmd = []

        if aCaseName != None:
            #find user case
            aCase = None
            aCases =  CFDSTUDYGUI_DataModel.GetCaseList(aStudy)
            for c in aCases:
                if c.GetName() == aCaseName:
                    aCase = c
                    break

            if aCase == None:
                if Trace(): print "Error. Can't find selected by user case"
                return

            # object of DATA folder
            aChildList = CFDSTUDYGUI_DataModel.ScanChildren(aCase, "^DATA$")
            if not len(aChildList) == 1:
                # no DATA folder
                if Trace(): print "There are not data folder in selected by user case"
                return

            aDataObj =  aChildList[0]
            aDataPath = CFDSTUDYGUI_DataModel._GetPath(aDataObj)
            if Trace(): print "Data Path: ", aDataPath

            # object of 'CFDSTUDYGUI' file
            if CFD_Code() == CFD_Saturne:
                aChildList = CFDSTUDYGUI_DataModel.ScanChildren(aDataObj, "^SaturneGUI$")
            elif CFD_Code() == CFD_Neptune:
                aChildList = CFDSTUDYGUI_DataModel.ScanChildren(aDataObj, "^NeptuneGUI$")
            if len(aChildList) == 0:
                # no 'CFDSTUDYGUI' file
                if Trace(): print "no 'SaturneGUI' or 'NeptuneGUI' file"
                return

        if CFD_Code() == CFD_Saturne:
            dlg = self.DialogCollector.GUIActivationDialog
            if dlg.ifUseLangOption():
                mess = "Language option not yet implemented "
                if QMessageBox.warning(None, "Warning", mess, QMessageBox.Yes, QMessageBox.No) == QMessageBox.No:
                    return

        if self.DialogCollector.GUIActivationDialog.isUseXmlFile() and \
               self.DialogCollector.GUIActivationDialog.currentXMLfile() != None:
            aXmlFile = str(self.DialogCollector.GUIActivationDialog.currentXMLfile())
            #XMLfile is chosen in the list into the dialog box
            if CFDSTUDYGUI_SolverGUI.findDockWindow(str(aXmlFile), aCase.GetName(),aStudy.GetName()) :
                self.updateActions()
                mess = aStudy.GetName() + " " + aCase.GetName() + " : " + aXmlFile + " is already launched"
                QMessageBox.information(None, "Information", mess, QMessageBox.Ok, QMessageBox.NoButton)
                return
            aCmd.append('-f')
            if CFD_Code() == CFD_Saturne:
                aCmd.append(aXmlFile)
            elif CFD_Code() == CFD_Neptune:
                import os.path
                aCmd.append(os.path.join("DATA", aXmlFile))

        else:
            aCmd.append('-n')

        if Trace(): print "CFDSTUDYGUI_ActionHandler:::::::: Args Command for ExecGui: ", aCmd
        wm = self._SolverGUI.ExecGUI(self.dskAgent().workspace(), aXmlFile, aCase, aCmd)

        self.updateActions()


    def __compile(self, aCaseObject):
        """
        Private method.
        Build the 'code_saturne compile -t' or the 'neptune_cfd compile -t' command.

        @type theCase: C{SObject}
        @param theCase: object from the Object Browser.
        @rtype: C{String}
        @return: command line
        """
        # object of SRC folder
        aChildList = CFDSTUDYGUI_DataModel.ScanChildren(aCaseObject, "SRC")
        if not len(aChildList) == 1:
            raise ValueError, "There is a mistake with the SRC directory"

        b, c = BinCode()
        cmd = b + " compile -t"

        return cmd


    def slotRunCase(self, study=None, case=None):
        """
        Run a case of Code_Saturne/NEPTUNE_CFD.

        @type study: C{String}
        @param study: name of the study.
        @type case: C{String}
        @param case: name of the case.
        """
        if Trace(): print 'slotRunCase'

        #select current case
        #get current selection
        sobj = self._singleSelectedObject()

        if Trace(): print "sobj: ", sobj

        if study:
            aStudy = study
        else:
            aStudy = CFDSTUDYGUI_DataModel.GetStudyByObj(sobj)
        if aStudy == None:
            if Trace(): print "Can't find CFD Study object"
            return

        # get current case
        if case:
            aCase = case
        else:
            aCase = CFDSTUDYGUI_DataModel.GetCase(sobj)
        if Trace(): print "CFDSTUDYGUI_ActionsHandler::::::::::::::::::: slotRunGUI Case: ", aCase
        self.DialogCollector.RunCaseDialog.setCurrentCase(aCase)
        self.DialogCollector.RunCaseDialog.setCurrentStudy(aStudy)
        self.DialogCollector.RunCaseDialog.fillData()
        self.DialogCollector.RunCaseDialog.exec_()
        if not self.DialogCollector.RunCaseDialog.result() == QDialog.Accepted:
            return

        print 'self.DialogCollector.RunCaseDialog.CompileModeBtn.isChecked()', self.DialogCollector.RunCaseDialog.CompileModeBtn.isChecked()

        if self.DialogCollector.RunCaseDialog.CompileModeBtn.isChecked():
            # get current case
            aCase = None
            aCaseName = str(self.DialogCollector.RunCaseDialog.CaseCB.currentText().toLatin1())
            aCaseList = CFDSTUDYGUI_DataModel.GetCaseList(aStudy)
            for c in aCaseList:
                if c.GetName() == aCaseName:
                    aCase = c
                    break

            if Trace(): print "Case: ", aCase
            if aCase == None:
                return

#            #self._CommandMgr.runFunctionDlg(self.__compile,
#            self._CommandMgr.runCommandDlg(self.__compile,
#                                           self.tr("STMSG_CHECK_COMPILATION"),\
#                                           False,\
#                                           aCaseObject = aCase)
            cmd = self.__compile(aCase)
            aChildList = CFDSTUDYGUI_DataModel.ScanChildren(aCase, "SRC")
            aSRCObj =  aChildList[0]
            aSRCPath = CFDSTUDYGUI_DataModel._GetPath(aSRCObj)

            dlg = CFDSTUDYGUI_CommandMgr.CFDSTUDYGUI_QProcessDialog(sgPyQt.getDesktop(),
                                                                    self.tr("STMSG_CHECK_COMPILATION"),
                                                                    [cmd],
                                                                    aSRCPath)
            dlg.show()
            if Trace(): print 'Compile done'
        else:
            #RunCase: activation of ics
            if Trace(): print 'RunMode'
            aCaseName = str(self.DialogCollector.RunCaseDialog.CaseCB.currentText().toLatin1())

            aCaseList = CFDSTUDYGUI_DataModel.GetCaseList(aStudy)
            for c in aCaseList:
                if c.GetName() == aCaseName:
                    aCase = c
                    break
            self.slotRunGUI(study=aStudy, case=aCase)


    def slotMeshConvertToMed(self):
        """
        """
        study = CFDSTUDYGUI_DataModel._getStudy()

        sg = CFDSTUDYGUI_DataModel.sg
        if sg.SelectedCount() <= 0:
            # no selection
            return
        if sg.SelectedCount() == 1:
            sobj = self._singleSelectedObject()
            medFile = str(QFileInfo(sobj.GetName()).baseName().toLatin1())
            self.DialogCollector.ECSConversionDialog.setResultFileName(medFile)
        else:
            self.DialogCollector.ECSConversionDialog.setResultFileName('')

        self.DialogCollector.ECSConversionDialog.exec_()
        if not self.DialogCollector.ECSConversionDialog.result() == QDialog.Accepted:
            return

        #curd = os.path.abspath('.')
        aFirtsObj = None
        if sg.SelectedCount() == 1:
            aFirtsObj = self._singleSelectedObject()
        else:
            list_obj  = self._multipleSelectedObject()
            if not list_obj == []:
                aFirtsObj = list_obj[0] #sg.getSelected(0)
            else:
                return

        aStudyObj = CFDSTUDYGUI_DataModel.GetStudyByObj(aFirtsObj)
        aChList = CFDSTUDYGUI_DataModel.ScanChildren(aStudyObj, "MESH")
        if not len(aChList) == 1:
            mess = "Directory MESH does not exist !"
            QMessageBox.critical(self, "Error", mess, QMessageBox.Ok, 0)
            return

        aMeshFold = aChList[0]
        thePath = CFDSTUDYGUI_DataModel._GetPath(aMeshFold)

        log.debug("slotMeshConvertToMed -> thePath = %s" % thePath)
        args = []

        b, c = BinCode()
        args.append(c)
        args.append("--mesh")
        if sg.SelectedCount() == 1:
            #args.append(sobj.GetName())
            args.append(CFDSTUDYGUI_DataModel._GetPath(sobj))
        else:
            for sobject in list_obj:
                args.append(CFDSTUDYGUI_DataModel._GetPath(sobject))

        outfile = self.DialogCollector.ECSConversionDialog.resultFileName()

        args.append("--sim-comm")
        args.append("--case")
        args.append(os.path.join(thePath, outfile))
        args.append("--med")
        if self.DialogCollector.ECSConversionDialog.IsVolOption():
            args.append("--volume")
        if self.DialogCollector.ECSConversionDialog.IsBordOption():
            args.append("--boundary")
        if self.DialogCollector.ECSConversionDialog.IsColor2GroupOption():
            args.append("--color-to-group")

        log.debug("slotMeshConvertToMed -> args = %s" % args)
        self._CommandMgr.runCommandDlg(aMeshFold,self.tr("STMSG_ECS_CONVERT"), args, thePath)

        self.updateObjBrowser(aMeshFold)


    def slotCopyCaseFile(self):
        """
        Copy data xml file from a study case to another with popup menu attached to the xml file
        """
        sobj = self._singleSelectedObject()
        if Trace(): print "sobj: ", sobj
        if sobj == None:
            return

        self.DialogCollector.CopyDialog.setCurrentObject(sobj)
        self.DialogCollector.CopyDialog.show()

        if not self.DialogCollector.CopyDialog.result() == QDialog.Accepted:
            return

        #update Object Browser
        aStudy = CFDSTUDYGUI_DataModel.GetStudyByObj(sobj)
        aCaseName = self.DialogCollector.CopyDialog.destCaseName()
        aCaseList = CFDSTUDYGUI_DataModel.GetCaseList(aStudy)
        for case in aCaseList:
            if aCaseName == case.GetName():
                self.updateObjBrowser(case)
                break

        # self.slotUpdateObjectBrowser()


    def slotCheckCompilation(self):
        """
        """
        #code for check of fortran files
        if Trace(): print 'slotCheckCompilation'

        #get current selection
        sobj = self._singleSelectedObject()
        if Trace(): print "sobj: ", sobj
        if sobj == None:
            return
        # get current case
        aCase = CFDSTUDYGUI_DataModel.GetCase(sobj)
        if Trace(): print "Case: ", aCase
        if aCase == None:
            return

#        self._CommandMgr.runFunctionDlg(self.__compile,
#                                        self.tr("STMSG_CHECK_COMPILATION"),\
#                                        False,\
#                                        aCaseObject = aCase)

        cmd = self.__compile(aCase)
        aChildList = CFDSTUDYGUI_DataModel.ScanChildren(aCase, "SRC")
        aSRCObj =  aChildList[0]
        aSRCPath = CFDSTUDYGUI_DataModel._GetPath(aSRCObj)

        dlg = CFDSTUDYGUI_CommandMgr.CFDSTUDYGUI_QProcessDialog(sgPyQt.getDesktop(),
                                                                self.tr("STMSG_CHECK_COMPILATION"),
                                                                [cmd],
                                                                aSRCPath)
        dlg.show()


    def slotRunScript(self):
        """
        """
        sobj = self._singleSelectedObject()
        if sobj:
            curd = os.path.abspath('.')
            father = sobj.GetFather()
            fatherpath = CFDSTUDYGUI_DataModel._GetPath(father)
            path = CFDSTUDYGUI_DataModel._GetPath(sobj)

            # check exec rights
            if not os.access(path, os.F_OK or os.X_OK):
                mess = self.tr("RUN_SCRIPT_ACTION_ACCESS_ERROR")
                QMessageBox.critical(None, "Error", mess, QMessageBox.Ok, 0)
                return

            father_father = father.GetFather()
            aChList = CFDSTUDYGUI_DataModel.ScanChildren(father_father, "RESU")
            if len(aChList) != 0:
                for i in aChList:
                    thePath = CFDSTUDYGUI_DataModel._GetPath(i)
                    if Trace(): print "CFDSTUDYGUI_ActionsHandler.slotRunScript::::::THE PATH", thePath
            self._CommandMgr.runCommandDlg(aChList[0], self.tr("STMSG_RUN_SCRIPT"), path, fatherpath)
            self.updateObjBrowser(aChList[0])


    def slotSaveDataFile(self):
        """
        Redirects B{Save} method to GUI of current solver
        """
        if Trace(): print "Call of GUI 'save' action"
        self._SolverGUI.onSaveXmlFile()


    def slotSaveAsDataFile(self):
        """
        Redirects B{SaveAs} method to GUI of current solver
        """
        if Trace(): print "Call of GUI 'save as' action"
        self._SolverGUI.onSaveAsXmlFile()


    def slotOpenShell(self):
        """
        Redirects B{OpenShell} method to GUI of current solver
        """
        self._SolverGUI.onOpenShell()


    def slotDisplayCurrentCase(self):
        """
        Redirects B{Display Current Case} method to GUI of current solver
        """
        self._SolverGUI.onDisplayCase()


    def slotHelpAbout(self):
        """
        Redirects B{About QDialog} display to GUI of current solver
        """
        self._SolverGUI.onHelpAbout()


    def slotSaturneReloadModule(self):
        """
        Redirects OpenShell method to GUI of current solver
        """
        self._SolverGUI.onSaturneReloadModule()


    def slotSaturneReloadPage(self):
        """
        Redirects OpenShell method to GUI of current solver
        """
        self._SolverGUI.onSaturneReloadPage()


    def slotSaturneHelpLicense(self):
        """
        Redirects OpenShell method to GUI of current solver
        """
        self._SolverGUI.onSaturneHelpLicense()


    def slotSaturneHelpCS(self):
        """
        Redirects OpenShell method to GUI of current solver
        """
        self._SolverGUI.onSaturneHelpCS()


    def slotSaturneHelpSD(self):
        """
        Redirects OpenShell method to GUI of current solver
        """
        self._SolverGUI.onSaturneHelpSD()


    def slotSaturneHelpCS_Kernel(self):
        """
        Redirects OpenShell method to GUI of current solver
        """
        self._SolverGUI.onSaturneHelpCS_Kernel()


    def slotSaturneHelpCS_Infos(self):
        """
        Redirects OpenShell method to GUI of current solver
        """
        self._SolverGUI.onSaturneHelpCS_Infos()


    def slotNeptuneWinBrowser(self, flag):
        """
        Redirects OpenShell method to GUI of current solver
        """
        self._SolverGUI.onNeptuneWinBrowser(flag)


    def slotNeptuneWinIdenty(self, flag):
        """
        Redirects OpenShell method to GUI of current solver
        """
        self._SolverGUI.onNeptuneWinIdenty(flag)


    def slotStopSolver(self):
        """
        Stops current solver process
        """
        sobj = self._singleSelectedObject()
        if sobj == None:
            return


    def slotShowSolverProcess(self):
        """
        Show current solver process
        """
        sobj = self._singleSelectedObject()
        if sobj == None:
            return


    def commonAction(self, theId):
        """
        Returns action by id from common action map of module
        """
        if not theId in self._CommonActionIdMap:
            raise ActionError, "Invalid action id"

        action_id = self._CommonActionIdMap[theId]

        if action_id == None or not action_id in self._ActionMap:
            raise ActionError, "Invalid action map content"
        return self._ActionMap[action_id]


    def solverAction(self, theId):
        """
        Returns action by id from solver action maps of module
        """
        action_id = None

        if theId in self._SolverActionIdMap:
            action_id =  self._SolverActionIdMap[theId]
        elif theId in self._SaturneActionIdMap:
            action_id = self._SaturneActionIdMap[theId]
        elif theId in self._NeptuneActionIdMap:
            action_id = self._NeptuneActionIdMap[theId]

        if action_id == None:
            raise ActionError, "Invalid action id"

        if not action_id in self._ActionMap:
            raise ActionError, "Invalid action map content"

        return self._ActionMap[action_id]


    def actionId(self, theId):
        """
        """
        action_id = None

        if theId in self._CommonActionIdMap:
            action_id =  self._CommonActionIdMap[theId]
        elif theId in self._SolverActionIdMap:
            action_id =  self._SolverActionIdMap[theId]
        elif theId in self._SaturneActionIdMap:
            action_id = self._SaturneActionIdMap[theId]
        elif theId in self._NeptuneActionIdMap:
            action_id = self._NeptuneActionIdMap[theId]

        if action_id == None:
            raise ActionError, "Invalid action id"

        return action_id


    def onCFDCode(self):
        """
        Shows/hides actual actions for current solver
        """
        import SalomePyQt
        if Trace(): print 'CFDSTUDYGUI_ActionsHandler.onCFDCode'

        is_cs = CFD_Code() == CFD_Saturne
        is_nc = CFD_Code() == CFD_Neptune

        #activate/deactivate saturne actions
        for a in self._SaturneActionIdMap:
            if a != SaturneHelpUserManualMenu:
                self.solverAction(a).setVisible(is_cs)

        #activate/deactivate neptune actions
        for a in self._NeptuneActionIdMap:
            if a != NeptuneWinMenu:
                self.solverAction(a).setVisible(is_nc)

        #update labels on common solver actions

        menuFileName = str(ObjectTR.tr("MEN_DESK_FILE"))
        menuHelpName = str(ObjectTR.tr("MEN_DESK_HELP"))

        menuFile = sgPyQt.getPopupMenu(menuFileName) #Info: return a QMenu
        menuHelp = sgPyQt.getPopupMenu(menuHelpName)

        menulisteFile_actions = menuFile.actions()
        menulisteHelp_actions = menuHelp.actions()

        if is_cs:
            codeFileName = str(ObjectTR.tr("SATURNE_FILE_MENU_TEXT"))
            codeHelpName = str(ObjectTR.tr("SATURNE_HELP_MENU_TEXT"))
            #sgPyQt.getPopupMenu(SalomePyQt.Window).setVisible(self.actionId(NeptuneWinMenu), False)
            #sgPyQt.getPopupMenu(SalomePyQt.Help).setVisible(self.actionId(SaturneHelpUserManualMenu), True)
        elif is_nc:
            codeFileName = str(ObjectTR.tr("NEPTUNE_FILE_MENU_TEXT"))
            codeHelpName = str(ObjectTR.tr("NEPTUNE_HELP_MENU_TEXT"))

        for i in range(len(menulisteFile_actions)):
            if menulisteFile_actions[i].text() == str(ObjectTR.tr("SOLVER_FILE_MENU_TEXT")):
                menulisteFile_actions[i].setText(codeFileName)
        for i in range(len(menulisteHelp_actions)):
            if menulisteHelp_actions[i].text() == str(ObjectTR.tr("SOLVER_HELP_MENU_TEXT")):
                menulisteHelp_actions[i].setText(codeHelpName)
                menulisteHelp_actions[i].setIconVisibleInMenu(True)

#    def processMgr(self):
#        """
#        Returns the process manager.
#        """
#        pass


    def dskAgent(self):
        """
        Returns the dekstop Agent.
        """
        return self._DskAgent


    def disconnectSolverGUI(self) :
        """
        Hide all the dock windows of CFDSTUDY, when activating another Salome Component
        We can have one or several of them with the right click on the main menu bar of
        Salome
        """
        #self._SolverGUI.disconnectDockWindowsBrowser()
        self._SolverGUI.disconnectDockWindows()



    def connectSolverGUI(self) :
        """
        Show all the dock windows of CFDSTUDY, when activating another Salome Component
        """
        self._SolverGUI.connectDockWindows()
        #self._SolverGUI.connectDockWindowsBrowser()
