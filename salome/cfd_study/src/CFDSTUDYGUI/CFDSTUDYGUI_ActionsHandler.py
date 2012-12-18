# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2012 EDF S.A.
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

import os, string, shutil
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

import SALOMEDS #Si on veut changer de couleur...
import salome,smesh
import SMESH

#-------------------------------------------------------------------------------
# Application modules
#-------------------------------------------------------------------------------

import CFDSTUDYGUI_DialogCollector
import CFDSTUDYGUI_DataModel
import CFDSTUDYGUI_Commons
import CFDSTUDYGUI_CommandMgr
from CFDSTUDYGUI_Agents import *
from CFDSTUDYGUI_Commons import CFD_Code, BinCode, CFD_Saturne, CFD_Neptune, sgPyQt, sg, CheckCFD_CodeEnv
import CFDSTUDYGUI_SolverGUI

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("CFDSTUDYGUI_ActionsHandler")
log.setLevel(logging.NOTSET)

#-------------------------------------------------------------------------------
# Global definitions
#-------------------------------------------------------------------------------

# Actions
SetStudyAction                = 1
AddCaseAction                 = 2
RunCaseAction                 = 3
LaunchGUIAction               = 4
OpenGUIAction                 = 5
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

#export/convert actions
ExportInPostProAction         = 40
ExportInSMESHAction           = 41
ConvertMeshToMed              = 42

#other actions
CheckCompilationAction        = 50
RunScriptAction               = 51

#Display Actions
DisplayMESHAction              = 60
DisplayGroupMESHAction         = 61
DisplayOnlyGroupMESHAction     = 62
HideGroupMESHAction            = 63
HideMESHAction                 = 64

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
SolverSaveAction               = 101
SolverSaveAsAction             = 102
SolverCloseAction              = 103
SolverUndoAction               = 104
SolverRedoAction               = 105

SolverToolsMenu                = 110
SolverOpenShellAction          = 111
SolverDisplayCurrentCaseAction = 112

SolverHelpMenu                 = 130
SolverHelpAboutAction          = 131 

#Help menu
SolverHelpLicense              = 251
SolverHelpGuidesMenu           = 260
SolverHelpUserGuide            = 261
SolverHelpTutorial             = 262
SolverHelpTheory               = 263
SolverHelpRefcard              = 264

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

        self.l_color = [(1,0,0),(0,1,0),(0,0,1),(1,1,0),(1,0,1),(0,1,1),]#(0.5,0,0),(0,0.5,0),(0,0,0.5),(0.2,0,0),(0,0.2,0),(0,0,0.2)]
        self.ul_color = []
        #intialise all dialogs
        self.DialogCollector = CFDSTUDYGUI_DialogCollector.CFDSTUDYGUI_DialogCollector()

        self._ActionMap = {}
        self._CommonActionIdMap = {}
        self._SolverActionIdMap = {}
        self._HelpActionIdMap = {}

        self._SalomeSelection = sgPyQt.getSelection()
        self._SolverGUI = CFDSTUDYGUI_SolverGUI.CFDSTUDYGUI_SolverGUI()
        self._DskAgent = Desktop_Agent()

        self.myVisu = None
        self.myViewManager = None

        try:
            import VISU
            import visu_gui
            self.myVisu = visu_gui.myVisu
            self.myViewManager = self.myVisu.GetViewManager()
            #self.myView = myViewManager.Create3DView()
            #self.myView = myViewManager.GetCurrentView()
            log.debug("__init__ myVisu = %s" % self.myVisu)
            log.debug("__init__ myViewManager = %s" % self.myViewManager)
            #log.debug("__init__ myView = %s" % self.myView)
        except Exception:
            log.debug("VISU module not available.")
            pass


    def createActions(self):
        """
        Creates menu, actions, and separators.
        """
        log.debug("createActions")
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
        #sgPyQt.createMenu(action, menu_id)
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
        # popup launch GUI on CFD CASE with slotRunGUI
        #sgPyQt.createMenu(action, menu_id)
        sgPyQt.createTool(action, tool_id)
        action_id = sgPyQt.actionId(action)
        self._ActionMap[action_id] = action
        self._CommonActionIdMap[LaunchGUIAction] = action_id
        self.connect(action, SIGNAL("activated()"), self.slotRunGUI)

        action = sgPyQt.createAction(-1,\
                                      ObjectTR.tr("Open GUI"),\
                                      ObjectTR.tr("LAUNCH_CFDSTUDY_GUI_TIP"),\
                                      ObjectTR.tr("LAUNCH_CFDSTUDY_GUI_SB"),\
                                      ObjectTR.tr("LAUNCH_CFDSTUDY_GUI_ICON"))
        action_id = sgPyQt.actionId(action)
        self._ActionMap[action_id] = action
        self._CommonActionIdMap[OpenGUIAction] = action_id
        self.connect(action, SIGNAL("activated()"), self.slotOpenCFD_GUI)

        action = sgPyQt.createAction(-1,\
                                      ObjectTR.tr("LAUNCH_CFDSTUDY_CASE_TEXT"),\
                                      ObjectTR.tr("LAUNCH_CFDSTUDY_CASE_TIP"),\
                                      ObjectTR.tr("LAUNCH_CFDSTUDY_CASE_SB"),\
                                      ObjectTR.tr("LAUNCH_CFDSTUDY_CASE_ICON"))
        # Run Case popup on study CFD name in object Browser with slotRunCase
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
        #sgPyQt.createTool(action, tool_id)
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

        if (self.myVisu != None and self.myViewManager != None):
            action = sgPyQt.createAction(-1,
                                          ObjectTR.tr("EXPORT_IN_POSTPRO_ACTION_TEXT"),
                                          ObjectTR.tr("EXPORT_IN_POSTPRO_ACTION_TIP"),
                                          ObjectTR.tr("EXPORT_IN_POSTPRO_ACTION_SB"),
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

# popup added to hide the mesh.

        action = sgPyQt.createAction(-1,\
                                      "Hide mesh",\
                                      "Hide mesh",\
                                      "Hide mesh",\
                                      ObjectTR.tr("MESH_OBJ_ICON"))
        self.connect(action, SIGNAL("activated()"), self.slotHideMESH)
        action_id = sgPyQt.actionId(action)
        self._ActionMap[action_id] = action
        self._CommonActionIdMap[HideMESHAction] = action_id

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
                                      ObjectTR.tr("ECS_CONVERT_ACTION_TEXT"),\
                                      ObjectTR.tr("ECS_CONVERT_ACTION_TIP"),\
                                      ObjectTR.tr("ECS_CONVERT_ACTION_SB"),\
                                      ObjectTR.tr("ECS_CONVERT_ACTION_ICON"))
        self.connect(action, SIGNAL("activated()"), self.slotMeshConvertToMed)
        action_id = sgPyQt.actionId(action)
        self._ActionMap[action_id] = action
        self._CommonActionIdMap[ConvertMeshToMed] = action_id

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
        # find the menu File into the Main Menu Bar of Salome
        fileId = sgPyQt.createMenu( ObjectTR.tr("MEN_DESK_FILE"), -1, -1)

        # create my menu into  menu File at position 7
        action_id = sgPyQt.createMenu(ObjectTR.tr("SOLVER_FILE_MENU_TEXT"), fileId, -1, 7, 1)
        self._SolverActionIdMap[SolverFileMenu] = action_id

        # warning: a Separator is a QMenu item (a trait)
        # create a separator after my menu in position 8
        action = sgPyQt.createSeparator()
        sgPyQt.createMenu(action, fileId, -1, 8, 1)

        # Save action
        action = sgPyQt.createAction(SolverSaveAction,\
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
        self._SolverActionIdMap[SolverSaveAction] = action_id

        # Save As action
        action = sgPyQt.createAction(-1,\
                                      ObjectTR.tr("SOLVER_SAVEAS_ACTION_TEXT"),\
                                      ObjectTR.tr("SOLVER_SAVEAS_ACTION_TIP"),\
                                      ObjectTR.tr("SOLVER_SAVEAS_ACTION_SB"),\
                                      ObjectTR.tr("SOLVER_SAVEAS_ACTION_ICON"),
                                      Qt.SHIFT+Qt.CTRL+Qt.Key_A)
        sgPyQt.createTool(action, tool_id)
        sgPyQt.createMenu(action, self._SolverActionIdMap[SolverFileMenu], 100)
        self.connect(action, SIGNAL("activated()"), self.slotSaveAsDataFile)

        action_id = sgPyQt.actionId(action)
        self._ActionMap[action_id] = action
        self._SolverActionIdMap[SolverSaveAsAction] = action_id
        action = sgPyQt.createSeparator()
        sgPyQt.createMenu(action, 1, 0, 2)

        # close GUI action
        action = sgPyQt.createAction(-1,\
                                      ObjectTR.tr("Close GUI"),\
                                      ObjectTR.tr("CLOSE_CFD_GUI_ACTION_TIP"),\
                                      ObjectTR.tr("CLOSE_CFD_GUI_ACTION_SB"),\
                                      ObjectTR.tr("CLOSE_CFD_GUI_ACTION_ICON"),
                                      Qt.SHIFT+Qt.CTRL+Qt.Key_W)
        sgPyQt.createTool(action, tool_id)
        action_id = sgPyQt.actionId(action)
        self._ActionMap[action_id] = action
        self._SolverActionIdMap[SolverCloseAction] = action_id
        self.connect(action, SIGNAL("activated()"), self.slotCloseCFD_GUI)

        # Add separator
        action = sgPyQt.createSeparator()
        sgPyQt.createTool(action, tool_id)

        # Undo action
        action = sgPyQt.createAction(-1, "Undo", "Undo", "Undo", \
                                      ObjectTR.tr("UNDO_CFD_GUI_ACTION_ICON"))
        sgPyQt.createTool(action, tool_id)
        action_id = sgPyQt.actionId(action)
        self._ActionMap[action_id] = action
        self._SolverActionIdMap[SolverUndoAction] = action_id
        self.connect(action, SIGNAL("activated()"), self.slotUndo)

        # Redo action
        action = sgPyQt.createAction(-1, "Redo", "Redo", "Redo", \
                                      ObjectTR.tr("REDO_CFD_GUI_ACTION_ICON"))
        sgPyQt.createTool(action, tool_id)
        action_id = sgPyQt.actionId(action)
        self._ActionMap[action_id] = action
        self._SolverActionIdMap[SolverRedoAction] = action_id
        self.connect(action, SIGNAL("activated()"), self.slotRedo)

        # Tools Menu
        action = sgPyQt.createSeparator()
        sgPyQt.createMenu(action, menu_id, 0, -1)

        action_id = sgPyQt.createMenu(ObjectTR.tr("SOLVER_TOOLS_MENU_TEXT"), menu_id)
        self._SolverActionIdMap[SolverToolsMenu] = action_id

        # Open shell action
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
        self._HelpActionIdMap[0] = action_id

        # Help menu: insert a Solver Menu Help to the Main Menu Help of Salome

        helpId = sgPyQt.createMenu( ObjectTR.tr("MEN_DESK_HELP"), -1, -1)
        #Info: Separator created at the end of the Menu Help (when we did not indicate a number)

        action = sgPyQt.createSeparator()
        sgPyQt.createMenu(action, helpId)
        #Info: Solver Help Menu created at the end of the Menu Help of Salome(when we did not indicate a number)
        action_id = sgPyQt.createMenu("Code_Saturne NEPTUNE_CFD", helpId)
        self._SolverActionIdMap[SolverHelpMenu] = action_id

        m = "About CFD"
        action = sgPyQt.createAction(-1, m, m, m)
        action_id = sgPyQt.actionId(action)
        self._ActionMap[action_id] = action
        self._SolverActionIdMap[SolverHelpAboutAction] = action_id
        sgPyQt.createMenu(action, self._SolverActionIdMap[SolverHelpMenu])
        self.connect(action, SIGNAL("activated()"), self.slotHelpAbout)
        self._ActionMap[action_id].setVisible(True)

        m = "License"
        action = sgPyQt.createAction(SolverHelpLicense, m, m, m)
        sgPyQt.createMenu(action, self._SolverActionIdMap[SolverHelpMenu])
        action_id = sgPyQt.actionId(action)
        self._ActionMap[action_id] = action
        self._HelpActionIdMap[SolverHelpLicense] = action_id
        self.connect(action, SIGNAL("activated()"), self.slotHelpLicense)

        # Guides menu
        action_id = sgPyQt.createMenu("Code_Saturne and NEPTUNE_CFD Guides", self._SolverActionIdMap[SolverHelpMenu])
        self._HelpActionIdMap[SolverHelpGuidesMenu] = action_id

        m = "User guide"
        action = sgPyQt.createAction(SolverHelpUserGuide, m, m, m)
        sgPyQt.createMenu(action, self._HelpActionIdMap[SolverHelpGuidesMenu])
        action_id = sgPyQt.actionId(action)
        self._ActionMap[action_id] = action
        self._HelpActionIdMap[SolverHelpUserGuide] = action_id
        self.connect(action, SIGNAL("activated()"), self.slotHelpUserGuide)

        m = "Tutorial"
        action = sgPyQt.createAction(SolverHelpTutorial, m, m, m)
        sgPyQt.createMenu(action, self._HelpActionIdMap[SolverHelpGuidesMenu])
        action_id = sgPyQt.actionId(action)
        self._ActionMap[action_id] = action
        self._HelpActionIdMap[SolverHelpTutorial] = action_id
        self.connect(action, SIGNAL("activated()"), self.slotHelpTutorial)

        m = "Theoretical guide"
        action = sgPyQt.createAction(SolverHelpTheory, m, m, m)
        sgPyQt.createMenu(action, self._HelpActionIdMap[SolverHelpGuidesMenu])
        action_id = sgPyQt.actionId(action)
        self._ActionMap[action_id] = action
        self._HelpActionIdMap[SolverHelpTheory] = action_id
        self.connect(action, SIGNAL("activated()"), self.slotHelpTheory)

        m = "Reference card"
        action = sgPyQt.createAction(SolverHelpRefcard, m, m, m)
        sgPyQt.createMenu(action, self._HelpActionIdMap[SolverHelpGuidesMenu])
        action_id = sgPyQt.actionId(action)
        self._ActionMap[action_id] = action
        self._HelpActionIdMap[SolverHelpRefcard] = action_id
        self.connect(action, SIGNAL("activated()"), self.slotHelpRefcard)

#        action_id = sgPyQt.createMenu(ObjectTR.tr("MESH_OR_GROUP_REPRESENTATION"), -1, -1)
#        self._CommonActionIdMap[SolverHelpGuidesMenu] = action_id

        #action = sgPyQt.createAction(-1,\
                                      #ObjectTR.tr("MESH_OR_GROUP_REPRESENTATION_SHADED"),\
                                      #ObjectTR.tr("MESH_OR_GROUP_REPRESENTATION_SHADED"),\
                                      #ObjectTR.tr("MESH_OR_GROUP_REPRESENTATION_SHADED"))
        #action_id = sgPyQt.actionId(action)
        #self._ActionMap[action_id] = action
        #self._CommonActionIdMap[DisplayTypeSHADED] = action_id
        #self.connect(action, SIGNAL("activated()"), self.slotDisplayTypeSHADED)

        #action = sgPyQt.createAction(-1,\
                                      #ObjectTR.tr("MESH_OR_GROUP_REPRESENTATION_WIREFRAME"),\
                                      #ObjectTR.tr("MESH_OR_GROUP_REPRESENTATION_WIREFRAME"),\
                                      #ObjectTR.tr("MESH_OR_GROUP_REPRESENTATION_WIREFRAME"))
        #action_id = sgPyQt.actionId(action)
        #self._ActionMap[action_id] = action
        #self._CommonActionIdMap[DisplayTypeWIREFRAME] = action_id
        #self.connect(action, SIGNAL("activated()"), self.slotDisplayTypeWIREFRAME)


    def updateActions(self):
        """
        Updates all action according with current selection and study states.
        This function connected to selection change signal.
        """
        log.debug("updateActions")
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
        if sobj != None:
            isStudy = CFDSTUDYGUI_DataModel.checkType(sobj, CFDSTUDYGUI_DataModel.dict_object["Study"])
            self.commonAction(AddCaseAction).setEnabled(isStudy)
            self.commonAction(RunCaseAction).setEnabled(isStudy)
            aStudy = CFDSTUDYGUI_DataModel.GetStudyByObj(sobj)
            aCase = CFDSTUDYGUI_DataModel.GetCase(sobj)

            if aStudy != None and aCase != None:
                self.commonAction(LaunchGUIAction).setEnabled(CFDSTUDYGUI_DataModel.checkCaseLaunchGUI(aCase))
                self.commonAction(OpenGUIAction).setEnabled(CFDSTUDYGUI_DataModel.checkCaseLaunchGUI(aCase))
            else:
                self.commonAction(LaunchGUIAction).setEnabled(False)

        #enable / disable solver actions
        isActivatedView = self._SolverGUI.isActive() # Main GUI Window is active

        for a in self._SolverActionIdMap:
            if a != SolverFileMenu and a != SolverToolsMenu and a != SolverHelpMenu:
                self.solverAction(a).setEnabled(isActivatedView)

        for a in self._HelpActionIdMap:
            if a != SolverHelpGuidesMenu:
                self.solverAction(a).setEnabled(isActivatedView)
                if CFD_Code() == CFD_Neptune:
                    self.solverAction(SolverHelpRefcard).setEnabled(False)

        if sobj != None:
            if CFDSTUDYGUI_DataModel.checkType(sobj, CFDSTUDYGUI_DataModel.dict_object["DATAfileXML"]):
                self.solverAction(SolverCloseAction).setEnabled(False)
                self.solverAction(SolverSaveAction).setEnabled(False)
                self.solverAction(SolverSaveAsAction).setEnabled(False)
                self.solverAction(SolverUndoAction).setEnabled(False)
                self.solverAction(SolverRedoAction).setEnabled(False)
                boo = True

                xmlName      = sobj.GetName()
                caseName     = sobj.GetFather().GetFather().GetName()
                studyCFDName = sobj.GetFather().GetFather().GetFather().GetName()
                dockName = self._SolverGUI.getDockTitleNameFromOB(studyCFDName,caseName,xmlName)

                listOfOpenSalomeStudies = CFDSTUDYGUI_DataModel._getlistOfOpenStudies()

                if listOfOpenSalomeStudies != [] and len(listOfOpenSalomeStudies) >= 1:
                    for nameSalomeStudy in CFDSTUDYGUI_DataModel._getlistOfOpenStudies():
                        studyId = CFDSTUDYGUI_DataModel._getStudy_Id(nameSalomeStudy)
                        if CFDSTUDYGUI_SolverGUI._c_CFDGUI.d_CfdCases != {}:
                            if CFDSTUDYGUI_SolverGUI._c_CFDGUI.d_CfdCases.has_key(studyId):
                                if CFDSTUDYGUI_SolverGUI._c_CFDGUI.d_CfdCases[studyId] != []:
                                    dockListe, dockListeWB = CFDSTUDYGUI_SolverGUI._c_CFDGUI.getDockListes(studyId)
                                    for dock in dockListe:
                                        if dockName == dock.windowTitle():
                                            self.commonAction(OpenGUIAction).setEnabled(False)
                                            if studyId != sgPyQt.getStudyId():
                                                self.solverAction(SolverCloseAction).setEnabled(False)
                                                self.solverAction(SolverSaveAction).setEnabled(False)
                                                self.solverAction(SolverSaveAsAction).setEnabled(False)
                                                self.solverAction(SolverUndoAction).setEnabled(False)
                                                self.solverAction(SolverRedoAction).setEnabled(False)
                                            else:
                                                self.solverAction(SolverCloseAction).setEnabled(True)
                                                self.solverAction(SolverSaveAction).setEnabled(True)
                                                self.solverAction(SolverSaveAsAction).setEnabled(True)
                                                self.solverAction(SolverUndoAction).setEnabled(True)
                                                self.solverAction(SolverRedoAction).setEnabled(True)
                                            boo = False
                if boo:
                    self.commonAction(OpenGUIAction).setEnabled(True)


    def customPopup(self, id, popup):
        """
        Callback for fill popup menu according current selection state.
        Function called by C{createPopupMenu} from CFDSTUDYGUI.py

        @type id: C{int}
        @param id: type of the branch tree slected in the Object Brower.
        @type popup: C{QPopupMenu}
        @param popup: popup menu from the Object Browser.
        """
        log.debug("customPopup")
        if id == CFDSTUDYGUI_DataModel.dict_object["Study"]:
            popup.addAction(self.commonAction(AddCaseAction))
            popup.addAction(self.commonAction(UpdateObjBrowserAction))
            popup.addAction(self.commonAction(RunCaseAction))
        elif id == CFDSTUDYGUI_DataModel.dict_object["Case"]:
            popup.addAction(self.commonAction(LaunchGUIAction))
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
            popup.addAction(self.commonAction(LaunchGUIAction))
        elif id == CFDSTUDYGUI_DataModel.dict_object["DATAfileXML"]:
            popup.addAction(self.commonAction(OpenGUIAction))
            popup.addAction(self.solverAction(SolverCloseAction))
            popup.addAction(self.commonAction(CopyCaseFileAction))
        elif id == CFDSTUDYGUI_DataModel.dict_object["SRCFolder"]:
            popup.addAction(self.commonAction(CheckCompilationAction))
        elif id == CFDSTUDYGUI_DataModel.dict_object["SRCFile"]:
            popup.addAction(self.commonAction(CheckCompilationAction))
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
        elif id == CFDSTUDYGUI_DataModel.dict_object["RESUSubFolder"]:
            popup.addAction(self.commonAction(RemoveAction))
        elif id == CFDSTUDYGUI_DataModel.dict_object["RESUSubErrFolder"]:
            popup.addAction(self.commonAction(RemoveAction))
        elif id == CFDSTUDYGUI_DataModel.dict_object["RESSRCFile"]:
            popup.addAction(self.commonAction(ViewAction))
        elif id == CFDSTUDYGUI_DataModel.dict_object["HISTFile"]:
            popup.addAction(self.commonAction(ViewAction))
            popup.addAction(self.commonAction(ExportInPostProAction))
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
        elif id == CFDSTUDYGUI_DataModel.dict_object["DESFile"] \
             or id == CFDSTUDYGUI_DataModel.dict_object["CGNSFile"] \
             or id == CFDSTUDYGUI_DataModel.dict_object["GeomFile"] \
             or id == CFDSTUDYGUI_DataModel.dict_object["CaseFile"] \
             or id == CFDSTUDYGUI_DataModel.dict_object["NeuFile"] \
             or id == CFDSTUDYGUI_DataModel.dict_object["MSHFile"] \
             or id == CFDSTUDYGUI_DataModel.dict_object["HexFile"] \
             or id == CFDSTUDYGUI_DataModel.dict_object["UnvFile"]:
            popup.addAction(self.commonAction(ConvertMeshToMed))
        elif id == CFDSTUDYGUI_DataModel.dict_object["MEDFile"]:
            popup.addAction(self.commonAction(ExportInSMESHAction))
        elif id == CFDSTUDYGUI_DataModel.dict_object["MESHFile"]:
            popup.addAction(self.commonAction(ViewAction))
        elif id == CFDSTUDYGUI_DataModel.dict_object["DATFile"]:
            popup.addAction(self.commonAction(EditAction))
        elif id == CFDSTUDYGUI_DataModel.dict_object["POSTFile"]:
            popup.addAction(self.commonAction(ViewAction))
        elif id == "VTKViewer":
            popup.addAction(self.commonAction(DisplayTypeSHADED))
            popup.addAction(self.commonAction(DisplayTypeWIREFRAME))
        else:

            for sobj in self._multipleSelectedObject():
                if sobj != None:
                    if sobj.GetFatherComponent().GetName() == "Mesh":
                        if sobj.GetFather().GetName() == "Mesh":
                            #Comment: mesh under Mesh module root in the Object browser

                            CFDSTUDYGUI_DataModel.SetAutoColor(sobj.GetFather())

                            for i in [DisplayMESHAction, HideMESHAction]:
                                popup.addAction(self.commonAction(i))
                                self.commonAction(i).setEnabled(True)

                        meshGroupObject, group = CFDSTUDYGUI_DataModel.getMeshFromGroup(sobj) # on teste et on recupere le groupe

                        if meshGroupObject <> None:
                            if len(self.l_color) == 0:
                                self.l_color = self.ul_color
                            if len(self.l_color) <> 0:
                                a = self.l_color[0]
                                self.ul_color.append(a)
                                self.l_color.remove(a)
                                x,y,z=a
                                group.SetColor(SALOMEDS.Color(x,y,z))

                            for i in [DisplayGroupMESHAction, DisplayOnlyGroupMESHAction, HideGroupMESHAction]:
                                popup.addAction(self.commonAction(i))
                                self.commonAction(i).setEnabled(True)


    def slotStudyLocation(self):
        """
        Loads the CFD study location. If the name of the CFD study
        does not exists, the corresponding folder is created.
        """
        log.debug("slotStudyLocation")
        dialog = self.DialogCollector.SetTreeLocationDialog
        dialog.exec_()
        if not self.DialogCollector.SetTreeLocationDialog.result() == QDialog.Accepted:
            return

        cursor = QCursor(Qt.BusyCursor)
        QApplication.setOverrideCursor(cursor)

        iok = CFDSTUDYGUI_DataModel._SetStudyLocation(theStudyPath = dialog.StudyPath,
                                                      theCaseNames = dialog.CaseNames)
        if iok:
            studyId = sgPyQt.getStudyId()
            sgPyQt.updateObjBrowser(studyId, 1)
            self.updateActions()

        QApplication.restoreOverrideCursor()


    def slotAddCase(self):
        """
        Builds new CFD cases.
        """
        log.debug("slotAddCase")
        dialog = self.DialogCollector.SetTreeLocationDialog
        dialog.setCaseMode(True)

        studyObj = self._singleSelectedObject()
        if studyObj == None:
            return

        dialog.StudyPath = CFDSTUDYGUI_DataModel._GetPath(studyObj)

        if not os.path.exists(dialog.StudyPath):
            mess = self.tr("ENV_DLG_INVALID_DIRECTORY").arg(dialog.StudyPath) + self.tr("STMSG_UPDATE_STUDY_INCOMING")
            QMessageBox.information(None, "Information", mess, QMessageBox.Ok, QMessageBox.NoButton)
            return

        dialog.exec_()
        if  not self.DialogCollector.SetTreeLocationDialog.result() == QDialog.Accepted:
            #cancel of new case creation
            dialog.setCaseMode(False)
            return

        dialog.setCaseMode(False)

        iok = CFDSTUDYGUI_DataModel._SetStudyLocation(theStudyPath = dialog.StudyPath,
                                                      theCaseNames = dialog.CaseNames)
        self.updateObjBrowser()


    def slotInfo(self):
        """
        Shows the QDialog with the info from CFDSTUDY:
            - CFD code selected
            - environnement variables defined
        """
        log.debug("slotInfo")
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
        log.debug("updateObjBrowser")
        cursor = QCursor(Qt.BusyCursor)
        QApplication.setOverrideCursor(cursor)

        CFDSTUDYGUI_DataModel.UpdateSubTree(Object)

        QApplication.restoreOverrideCursor()


    def slotViewAction(self):
        """
        Edits in the read only mode the file selected in the Object Browser.
        Warning, the editor is always emacs!
        """
        viewerName = str( sgPyQt.stringSetting( "CFDSTUDY", "ExternalEditor", self.tr("CFDSTUDY_PREF_EDITOR")).trimmed() )
        if viewerName != "":
            sobj = self._singleSelectedObject()
            if sobj is not None:
                path = CFDSTUDYGUI_DataModel._GetPath(sobj)
                if re.match(".*emacs$", viewerName):
                    os.spawnlp(os.P_NOWAIT, viewerName , viewerName, path, "-f", "toggle-read-only")
                elif re.match("vi", viewerName) or re.match("vim", viewerName):
                    os.system("xterm -sb -e vi "  + path )
                else:
                    os.spawnlp(os.P_NOWAIT, viewerName ,viewerName , path)


    def slotEditAction(self):
        """
        Edits in the user's editor the file selected in the Object Browser.
        """
        viewerName = str( sgPyQt.stringSetting( "CFDSTUDY", "ExternalEditor", self.tr("CFDSTUDY_PREF_EDITOR") ).trimmed() )

        if viewerName != "":
            #viewerName = str(viewer.toLatin1())
            sobj = self._singleSelectedObject()
            if not sobj == None:
                path = CFDSTUDYGUI_DataModel._GetPath(sobj)
                os.spawnlp(os.P_NOWAIT,viewerName ,viewerName , path)


    def slotRemoveAction(self):
        """
        Deletes file or folder from the Object Browser, and from the unix system files.
        Delete dock windows attached to a CFD Study if this study is deleted from the Object Browser.
        """
        log.debug("slotRemoveAction")
        sobj = self._singleSelectedObject()
        if not sobj == None:
            mess = ObjectTR.tr("REMOVE_ACTION_CONFIRM_MESS").arg(sobj.GetName())
            if QMessageBox.warning(None, "Warning", mess, QMessageBox.Yes, QMessageBox.No) == QMessageBox.No:
                return

            path = CFDSTUDYGUI_DataModel._GetPath(sobj)
            c = CFDSTUDYGUI_DataModel.GetCase(sobj).GetName()
            caseName  = sobj.GetName()
            studyName = CFDSTUDYGUI_DataModel.GetStudyByObj(sobj)

            father = sobj.GetFather()
            fatherpath = CFDSTUDYGUI_DataModel._GetPath(father)
            fathername = father.GetName()

            if c == caseName:
                try:
                    self._SolverGUI.removeDockWindow(studyName, caseName)
                except:
                    pass

            watchCursor = QCursor(Qt.WaitCursor)
            QApplication.setOverrideCursor(watchCursor)
            if os.path.isdir(path):
                shutil.rmtree(path)
            elif os.path.isfile(path):
                os.remove(path)
            QApplication.restoreOverrideCursor()

            self.updateObjBrowser(father)


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

            # chechpoint and mesh_input
            #if attr.Value() == CFDSTUDYGUI_DataModel.dict_object["PRETFolder"] or \
               #attr.Value() == CFDSTUDYGUI_DataModel.dict_object["SUITEFolder"]:
                #case = CFDSTUDYGUI_DataModel.GetCase(sobj)
                #if not case == None:
                    #iter  = study.NewChildIterator(case)
                    #while iter.More():
                        #if iter.Value().GetName() == "DATA":
                            #newpath = os.path.join(CFDSTUDYGUI_DataModel._GetPath(iter.Value()), sobj.GetName())

                            ##remove if exists
                            #if os.path.exists(newpath):
                                #mess = ObjectTR.tr("OVERWRITE_CONFIRM_MESS").arg(sobj.GetName())
                                #if QMessageBox.warning(None, "Warning", mess, QMessageBox.Yes, QMessageBox.No) == QMessageBox.No:
                                    #return
                                #os.remove(newpath)

                            #shutil.copy2(path, CFDSTUDYGUI_DataModel._GetPath(iter.Value()))
                            #CFDSTUDYGUI_DataModel._RebuildTreeRecursively(case)
                            #break
                        #iter.Next()
                    #sgPyQt.updateObjBrowser(studyId,1)
            #else:
            parent = sobj.GetFather()
            if not parent == None:
                parent = parent.GetFather()
                if not parent == None and parent.GetName() == "DATA":
                    parentPath = CFDSTUDYGUI_DataModel._GetPath(parent)
                    newpath = os.path.join(parentPath, sobj.GetName())

                    if os.path.exists(newpath):
                        mess = ObjectTR.tr("OVERWRITE_CONFIRM_MESS")
                        if QMessageBox.warning(None, "Warning", mess, QMessageBox.Yes, QMessageBox.No) == QMessageBox.No:
                            return

                    shutil.copy2(path, parent)
                    self.updateObjBrowser(parent)


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
                        destPath = os.path.join(parentPath, sobj.GetName())
                        if os.path.exists(destPath):
                            mess = ObjectTR.tr("OVERWRITE_CONFIRM_MESS")
                            if QMessageBox.warning(None, "Warning", mess, QMessageBox.Yes, QMessageBox.No) == QMessageBox.No:
                                return
                        shutil.copy2(path, parentPath)
                        self.updateObjBrowser(parent)


    def slotMoveToDRAFT(self):
        """
        """
        sobj = self._singleSelectedObject()
        if not sobj == None:
            path = CFDSTUDYGUI_DataModel._GetPath(sobj)
            parent = sobj.GetFather()
            if not parent == None:
                parentPath = os.path.join(CFDSTUDYGUI_DataModel._GetPath(parent), 'DRAFT')
                destPath = os.path.join(parentPath, sobj.GetName())
                if os.path.exists(destPath):
                    mess = ObjectTR.tr("OVERWRITE_CONFIRM_MESS")
                    if QMessageBox.warning(None, "Warning", mess, QMessageBox.Yes, QMessageBox.No) == QMessageBox.No:
                        return
                    else:
                        os.remove(destPath)

                if os.path.exists(parentPath) == False:
                    os.mkdir(parentPath)

                shutil.move(path, parentPath)
                self.updateObjBrowser(parent)


    def _singleSelectedObject(self):
        """
        """
        study = CFDSTUDYGUI_DataModel._getStudy()
        if sg.SelectedCount() == 1:
            entry = sg.getSelected(0)
            if entry != '':
                return study.FindObjectID(entry)
        return None


    def _multipleSelectedObject(self):
        """
        """
        study = CFDSTUDYGUI_DataModel._getStudy()

        i = 0
        liste_SObj = []
        while i < sg.SelectedCount():
            entry = sg.getSelected(i)
            if entry != '':
                liste_SObj.append(study.FindObjectID(entry))
            i = i+1
        return liste_SObj


    def slotExportInPostPro(self):
        """
        """
        waitCursor = QCursor(Qt.WaitCursor)
        QApplication.setOverrideCursor(waitCursor)

        sobj = self._singleSelectedObject()
        if not sobj == None:

            path = CFDSTUDYGUI_DataModel._GetPath(sobj)
            if re.match(".*\.med$", sobj.GetName()):
                #export Med file
                self.myVisu.ImportFile(path)
            elif re.match(".*\.dat$", sobj.GetName()) or re.match(".*\.csv$", sobj.GetName()):
                self.myVisu.ImportTables(path, True)
            studyId = sgPyQt.getStudyId()
            sgPyQt.updateObjBrowser(studyId,1)
        QApplication.restoreOverrideCursor()


    def slotExportInSMESH(self):
        """
        smesh_component         is a smeshDC.smeshDC instance
        SO_SMESH_COMPONENT         is a SALOMEDS._objref_SComponent instance
        aMeshes                 is a list of smeshDC.Mesh instances of the meshes into the med file
        meshDC.GetMesh()         is a Corba SMESH._objref_SMESH_Mesh instance
        SO_SMESH                 is a SALOMEDS._objref_SObject instance representing mesh object into
                                      Object browser under SMESH Component
        Create Med structure of the med file whose complete name is path,
             into smesh component and puplication of the mesh into Object Browser
             aMeshes is a list of smeshDC.Mesh instances of the meshes into the med file
             (we can have several meshes into a med file)

        """
        waitCursor = QCursor(Qt.WaitCursor)
        QApplication.setOverrideCursor(waitCursor)

        sobj = self._singleSelectedObject()
        if not sobj == None:
            path = CFDSTUDYGUI_DataModel._GetPath(sobj)
            studyId = salome.sg.getActiveStudyId()
            if smesh and re.match(".*\.med$", sobj.GetName()):
                smesh.SetCurrentStudy(salome.myStudy)
                aMeshes, aStatus = smesh.CreateMeshesFromMED(path)
                if not aStatus:
                    QApplication.restoreOverrideCursor()
                    mess = ObjectTR("EXPORT_IN_SMESH_ACTION_WARNING")
                    QMessageBox.warning(None, "Warning", mess, QMessageBox.Ok, 0)
                    return

                (reppath,fileName)=   os.path.split(path)
                for aMeshDC in aMeshes:
                    aMeshDC.SetAutoColor(1)
                    mesh = aMeshDC.GetMesh()

            sgPyQt.updateObjBrowser(studyId, 1)

        QApplication.restoreOverrideCursor()


    def slotDisplayMESH(self):
        """
        Changed on November 2010 for the popup menu: SMESH Mesh objects can have the slotDisplayMESH directly
        the old code with referenced objects is deleted
        """
        waitCursor = QCursor(Qt.WaitCursor)
        QApplication.setOverrideCursor(waitCursor)

        if self._multipleSelectedObject() == None:
            mess = "Display MESH: No object selected into Object Browser"
            QMessageBox.warning(None, "Warning", mess, QMessageBox.Ok, 0)
            return
        smeshgui = salome.ImportComponentGUI("SMESH")
        studyId = salome.sg.getActiveStudyId()
        smeshgui.Init(studyId)

        log.debug("slotDisplayMESH -> self._multipleSelectedObject()[0].GetName()= %s" % self._multipleSelectedObject()[0].GetName())
        for  sobj in self._multipleSelectedObject():
            if sobj != None:
                entry = sobj.GetID()
                if entry == None:
                    mess = "slotDisplayMESH: No mesh with the Name: " + sobj.GetName() + ", under Mesh into Object Browser"
                    QMessageBox.warning(None, "Warning", mess, QMessageBox.Ok, 0)
                    QApplication.restoreOverrideCursor()
                    return
                #Displaying Mesh
                if CFDSTUDYGUI_DataModel.getMeshFromMesh(sobj):
                    smeshgui.CreateAndDisplayActor(entry)
                    sgPyQt.updateObjBrowser(studyId,1)
                    salome.sg.UpdateView()
                    salome.sg.FitAll()
            else:
                mess = "slotDisplayMESH: Entry Id not stored for the mesh: " + sobj.GetName()
                QMessageBox.warning(None, "Warning", mess, QMessageBox.Ok, 0)
        QApplication.restoreOverrideCursor()


    def slotDisplayMESHGroups(self):
        """
        Changed on November 2010 for the popup menu: SMESH Group Mesh objects can have the slotDisplayMESHGroups directly
        the old code with referenced objects is deleted
        """
        waitCursor = QCursor(Qt.WaitCursor)
        QApplication.setOverrideCursor(waitCursor)

        if self._multipleSelectedObject() == None:
            mess = "Display MESH Groups: No object selected into Object Browser"
            QMessageBox.warning(None, "Warning", mess, QMessageBox.Ok, 0)
            return
        smeshgui = salome.ImportComponentGUI("SMESH")

        for sobj_group in self._multipleSelectedObject():
            if sobj_group != None:
                meshgroup,group = CFDSTUDYGUI_DataModel.getMeshFromGroup(sobj_group)
                if meshgroup:
                    smeshgui.CreateAndDisplayActor(sobj_group.GetID())
            else:
                mess = "No group "+ sobj_group.GetName() + " whose mesh father name is:",sobj_group.GetFatherComponent().GetName() #GetFather().GetFather().GetName()
                QMessageBox.warning(None, "Warning", mess, QMessageBox.Ok, 0)
        salome.sg.UpdateView()
        salome.sg.FitAll()

        QApplication.restoreOverrideCursor()


    def slotDisplayOnlyMESHGroups(self):
        """
        """
        waitCursor = QCursor(Qt.WaitCursor)
        QApplication.setOverrideCursor(waitCursor)

        sobj = self._singleSelectedObject()
        id = sobj.GetID()
        if id:
            salome.sg.EraseAll()
            #salome.sg.Display(entryIdGroup)#Only(entryIdGroup)
            smeshgui = salome.ImportComponentGUI("SMESH")
            smeshgui.CreateAndDisplayActor(id)
        else:
            mess = "No Entry Id for group "+ sobj.GetName() + " whose mes Name is:",sobj.GetFatherComponent().GetName() #GetFather().GetFather().GetName()
            QMessageBox.warning(None, "Warning", mess, QMessageBox.Ok, 0)
        salome.sg.UpdateView()
        salome.sg.FitAll()

        QApplication.restoreOverrideCursor()


    def slotHideMESHGroups(self):
        """
        """
        waitCursor = QCursor(Qt.WaitCursor)
        QApplication.setOverrideCursor(waitCursor)

        if self._multipleSelectedObject() == None:
            mess = "Hide MESH Groups: No object selected into Object Browser"
            QMessageBox.warning(None, "Warning", mess, QMessageBox.Ok, 0)
            return

        for sobj in self._multipleSelectedObject():
            id = sobj.GetID()
            if id:
                meshgroup,group = CFDSTUDYGUI_DataModel.getMeshFromGroup(sobj)
                if meshgroup:
                    salome.sg.Erase(id)
            else:
                mess = "No Entry Id for group "+ sobj.GetName() + " whose mesh Name is:",sobj.GetFatherComponent().GetName() # sobj.GetFather().GetFather().GetName()
                QMessageBox.warning(None, "Warning", mess, QMessageBox.Ok, 0)
        salome.sg.UpdateView()
        salome.sg.FitAll()

        QApplication.restoreOverrideCursor()


    def slotHideMESH(self):
        """
        Changed on November 2010 for the popup menu: SMESH Mesh objects can have the slotHideMESH directly
        """
        waitCursor = QCursor(Qt.WaitCursor)
        QApplication.setOverrideCursor(waitCursor)

        if self._multipleSelectedObject() == None:
            mess = "Hide MESH: No object selected into Object Browser"
            QMessageBox.warning(None, "Warning", mess, QMessageBox.Ok, 0)
            return

        for sobj in self._multipleSelectedObject():
            id = sobj.GetID()
            if id:
                if CFDSTUDYGUI_DataModel.getMeshFromMesh(sobj):
                    salome.sg.Erase(id)
            else:
                mess = "No Entry Id for mesh "+ sobj.GetName()
                QMessageBox.warning(None, "Warning", mess, QMessageBox.Ok, 0)
        salome.sg.UpdateView()
        salome.sg.FitAll()

        QApplication.restoreOverrideCursor()


    def OpenCFD_GUI(self,sobj):
        """
        Open into Salome the CFD GUI from an XML file whose name is sobj.GetName()
        """
        log.debug("OpenCFD_GUI")
        import os
        if sobj != None:
            if not os.path.exists(CFDSTUDYGUI_DataModel._GetPath(sobj)):
                mess = self.tr("ENV_DLG_INVALID_FILE").arg("CFD_Code").arg(CFDSTUDYGUI_DataModel._GetPath(sobj))+ self.tr("STMSG_UPDATE_STUDY_INCOMING")
                QMessageBox.information(None, "Information", mess, QMessageBox.Ok, QMessageBox.NoButton)
                self.updateObjBrowser()
                return
            if CFDSTUDYGUI_DataModel.checkType(sobj, CFDSTUDYGUI_DataModel.dict_object["DATAfileXML"]):
                aXmlFileName = sobj.GetName()
                aCase = CFDSTUDYGUI_DataModel.GetCase(sobj)
                aStudy = CFDSTUDYGUI_DataModel.GetStudyByObj(sobj)
                if aCase:
                    aCaseName = aCase.GetName()
                else:
                    mess = "Error: "+ aXmlFileName + " file has no CFD Case into the Salome Object browser"
                    QMessageBox.warning(None, "Warning", mess, QMessageBox.Ok, 0)
                    return
                if aStudy:
                    aStudyName = aStudy.GetName()
                else:
                    mess = "Error: "+ aXmlFileName + " file has no CFD Study into the Salome Object browser"
                    QMessageBox.warning(None, "Warning", mess, QMessageBox.Ok, 0)
                    return
                if CFDSTUDYGUI_SolverGUI.findDockWindow(aXmlFileName, aCaseName,aStudyName):
                    mess = aStudyName + " " + aCaseName + ": " + aXmlFileName + " is already opened"
                    QMessageBox.information(None, "Information", mess, QMessageBox.Ok, QMessageBox.NoButton)
                    return

                # xml case file not already opened
                aCmd = []
                aCmd.append('-p')
                aCmd.append(aXmlFileName)
                aXmlFile = sobj
                wm = self._SolverGUI.ExecGUI(self.dskAgent().workspace(), aXmlFile, aCase, aCmd)
                self.updateActions()


    def slotOpenCFD_GUI(self):
        """
        Open into Salome the CFD GUI from an XML file whose name is sobj.GetName()
        """
        log.debug("slotOpenCFD_GUI")
        sobj = self._singleSelectedObject()
        if sobj != None:
            import os
            if not os.path.exists(CFDSTUDYGUI_DataModel._GetPath(sobj)):
                mess = self.tr("ENV_DLG_INVALID_FILE").arg("CFD_Code").arg(CFDSTUDYGUI_DataModel._GetPath(sobj))+ self.tr("STMSG_UPDATE_STUDY_INCOMING")
                QMessageBox.information(None, "Information", mess, QMessageBox.Ok, QMessageBox.NoButton)
                self.updateObjBrowser()
                return
            self.OpenCFD_GUI(sobj)


    def CloseCFD_GUI(self, sobj):
        """
        Close into Salome the CFD GUI from an XML file whose name is sobj.GetName()
        """
        log.debug("CloseCFD_GUI")
        if sobj != None and CFDSTUDYGUI_DataModel.checkType(sobj, CFDSTUDYGUI_DataModel.dict_object["DATAfileXML"]):
            aXmlFileName = sobj.GetName()
            aCase = CFDSTUDYGUI_DataModel.GetCase(sobj)
            aStudy = CFDSTUDYGUI_DataModel.GetStudyByObj(sobj)
            if aCase:
                aCaseName = aCase.GetName()
            else:
                mess = "Error: "+ aXmlFileName + " file has no CFD Case into the Salome Object browser"
                QMessageBox.warning(None, "Warning", mess, QMessageBox.Ok, 0)
                return

            if aStudy:
                aStudyName = aStudy.GetName()
            else:
                mess = "Error: "+ aXmlFileName + " file has no CFD Study into the Salome Object browser"
                QMessageBox.warning(None, "Warning", mess, QMessageBox.Ok, 0)
                return
        else:
            # close the active CFDGUI window with the icon button CLOSE_CFD_GUI_ACTION_ICON in the tool bar
            aStudyName, aCaseName, aXmlFileName = self._SolverGUI.getStudyCaseXmlNames(self._SolverGUI._CurrentWindow)

        log.debug("CloseCFD_GUI %s %s %s" % (aStudyName, aCaseName, aXmlFileName))
        if self._SolverGUI.okToContinue():
            self._SolverGUI.removeDockWindow(aStudyName, aCaseName, aXmlFileName)
            self.commonAction(OpenGUIAction).setEnabled(True)
            self.updateActions()


    def slotCloseCFD_GUI(self):
        """
        Close into Salome the CFD GUI from an XML file whose name is sobj.GetName()
        """
        log.debug("slotCloseCFD_GUI")
        sobj = self._singleSelectedObject()
        self.CloseCFD_GUI(sobj)


    def slotUndo(self):
        self._SolverGUI.onUndo()


    def slotRedo(self):
        self._SolverGUI.onRedo()


    def slotRunGUI(self, study=None, case=None):
        """
        Build the command line for the GUI of Code_Saturne/NEPTUNE_CFD.

        @type study: C{String}
        @param study: name of the study.
        @type case: C{String}
        @param case: name of the case.
        """
        log.debug("slotRunGUI")
        #get current selection
        sobj = self._singleSelectedObject()
        import os
        if not os.path.exists(CFDSTUDYGUI_DataModel._GetPath(sobj)):
            mess = self.tr("ENV_DLG_INVALID_DIRECTORY").arg(CFDSTUDYGUI_DataModel._GetPath(sobj))+ self.tr("STMSG_UPDATE_STUDY_INCOMING")
            QMessageBox.information(None, "Information", mess, QMessageBox.Ok, QMessageBox.NoButton)
            self.updateObjBrowser()
            return
        if study:
            aStudy = study
        else:
            aStudy = CFDSTUDYGUI_DataModel.GetStudyByObj(sobj)

        if aStudy == None:
            return
        objname = sobj.GetName()

        # get current case
        if case:
            aCase = case
        else:
            aCase = CFDSTUDYGUI_DataModel.GetCase(sobj)

        self.DialogCollector.GUIActivationDialog.setCurrentCase(aCase)
        self.DialogCollector.GUIActivationDialog.setCurrentStudy(aStudy)

        if CFDSTUDYGUI_DataModel.checkType(sobj, CFDSTUDYGUI_DataModel.dict_object["DATAfileXML"]):
            xmlFileName = objname
        else:
            xmlFileName = ""

        self.DialogCollector.GUIActivationDialog.fillData(xmlFileName)

        if aCase != None:
            if CFDSTUDYGUI_DataModel.checkType(sobj, CFDSTUDYGUI_DataModel.dict_object["DATAfileXML"]) or \
               CFDSTUDYGUI_DataModel.checkType(sobj, CFDSTUDYGUI_DataModel.dict_object["RESXMLFile"]):
                #checks that such tab already opened
                # launch GUI from an xml file from the Object browser
                if CFDSTUDYGUI_SolverGUI.findDockWindow(sobj.GetName(), aCase.GetName(), aStudy.GetName()):
                    mess = aStudy.GetName() + " " + aCase.GetName() + ": " + sobj.GetName() + " is already opened"
                    #QMessageBox.information(None, "Information", mess, QMessageBox.Ok, QMessageBox.NoButton)
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
                return

            # object of DATA folder
            aChildList = CFDSTUDYGUI_DataModel.ScanChildren(aCase, "^DATA$")
            if not len(aChildList) == 1:
                # no DATA folder
                return

            aDataObj =  aChildList[0]
            aDataPath = CFDSTUDYGUI_DataModel._GetPath(aDataObj)
            # object of 'CFDSTUDYGUI' file
            if CFD_Code() == CFD_Saturne:
                aChildList = CFDSTUDYGUI_DataModel.ScanChildren(aDataObj, "^SaturneGUI$")
            elif CFD_Code() == CFD_Neptune:
                aChildList = CFDSTUDYGUI_DataModel.ScanChildren(aDataObj, "^NeptuneGUI$")
            if len(aChildList) == 0:
                # no 'CFDSTUDYGUI' file
                return

        xmlFiles = CFDSTUDYGUI_DataModel.ScanChildren(aDataObj, ".*")
        if self.DialogCollector.GUIActivationDialog.isUseXmlFile() and \
               self.DialogCollector.GUIActivationDialog.currentXMLfile() != None:
            aXmlFile = str(self.DialogCollector.GUIActivationDialog.currentXMLfile())
            #XMLfile is chosen in the list into the dialog box
            if CFDSTUDYGUI_SolverGUI.findDockWindow(str(aXmlFile), aCase.GetName(),aStudy.GetName()):
                self.updateActions()
                mess = aStudy.GetName() + " " + aCase.GetName() + ": " + aXmlFile + " is already launched"
                QMessageBox.information(None, "Information", mess, QMessageBox.Ok)
                return
            aCmd.append('-p')
            aCmd.append(aXmlFile)
        else:
            aCmd.append('-n')

        sobjxml = None
        if aXmlFile != None:
            for xmlf in xmlFiles:
                if xmlf.GetName() == aXmlFile:
                    sobjxml = xmlf
        if not os.path.exists(os.path.join(CFDSTUDYGUI_DataModel._GetPath(aCase),"DATA")):
            mess = self.tr("ENV_DLG_INVALID_DIRECTORY").arg(os.path.join(CFDSTUDYGUI_DataModel._GetPath(aCase),"DATA"))+ self.tr("STMSG_UPDATE_STUDY_INCOMING")
            QMessageBox.information(None, "Information", mess, QMessageBox.Ok, QMessageBox.NoButton)
            self.updateObjBrowser()
            return
        wm = self._SolverGUI.ExecGUI(self.dskAgent().workspace(), sobjxml, aCase, aCmd)
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
        cmd = ""
        aChildList = CFDSTUDYGUI_DataModel.ScanChildren(aCaseObject, "SRC")
        if not len(aChildList) == 1:
            raise ValueError, "There is a mistake with the SRC directory"

        env_code, mess = CheckCFD_CodeEnv(CFD_Code())

        if not env_code:
            QMessageBox.critical(self,"Error", mess, QMessageBox.Ok, 0)
        else:
            b, c,mess = BinCode()
            if mess == "":
                cmd = b + " compile -t"
            else:
                QMessageBox.critical(self,"Error", mess, QMessageBox.Ok, 0)
        return cmd


    def slotRunCase(self, study=None, case=None):
        """
        Run a case of Code_Saturne/NEPTUNE_CFD.

        @type study: C{String}
        @param study: name of the study.
        @type case: C{String}
        @param case: name of the case.
        """
        log.debug("slotRunCase")
        #select current case
        #get current selection
        sobj = self._singleSelectedObject()

        if study:
            aStudy = study
        else:
            aStudy = CFDSTUDYGUI_DataModel.GetStudyByObj(sobj)

        if aStudy == None:
            return

        # get current case
        if case:
            aCase = case
        else:
            aCase = CFDSTUDYGUI_DataModel.GetCase(sobj)

        self.DialogCollector.RunCaseDialog.setCurrentCase(aCase)
        self.DialogCollector.RunCaseDialog.setCurrentStudy(aStudy)
        self.DialogCollector.RunCaseDialog.fillData()
        self.DialogCollector.RunCaseDialog.exec_()

        if not self.DialogCollector.RunCaseDialog.result() == QDialog.Accepted:
            return

        if self.DialogCollector.RunCaseDialog.CompileModeBtn.isChecked():
            # get current case
            aCase = None
            aCaseName = str(self.DialogCollector.RunCaseDialog.CaseCB.currentText().toLatin1())
            aCaseList = CFDSTUDYGUI_DataModel.GetCaseList(aStudy)
            for c in aCaseList:
                if c.GetName() == aCaseName:
                    aCase = c
                    break

            if aCase == None:
                return

            cmd = self.__compile(aCase)
            if cmd != "":
                aChildList = CFDSTUDYGUI_DataModel.ScanChildren(aCase, "SRC")
                aSRCObj =  aChildList[0]
                aSRCPath = CFDSTUDYGUI_DataModel._GetPath(aSRCObj)

                dlg = CFDSTUDYGUI_CommandMgr.CFDSTUDYGUI_QProcessDialog(sgPyQt.getDesktop(),
                                                                        self.tr("STMSG_CHECK_COMPILATION"),
                                                                        [cmd],
                                                                        None,
                                                                        aSRCPath)
                dlg.show()
        else:
            #RunCase: activation of ics
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
        if sg.SelectedCount() != 1:
            # no selection
            return
        elif sg.SelectedCount() == 1:
            sobj = self._singleSelectedObject()
            medFile = str(QFileInfo(sobj.GetName()).baseName().toLatin1())
            self.DialogCollector.ECSConversionDialog.setResultFileName(medFile)

        self.DialogCollector.ECSConversionDialog.exec_()
        if not self.DialogCollector.ECSConversionDialog.result() == QDialog.Accepted:
            return

        aFirtsObj = None
        if sg.SelectedCount() == 1:
            aFirtsObj = self._singleSelectedObject()
        else:
            list_obj  = self._multipleSelectedObject()
            if not list_obj == []:
                aFirtsObj = list_obj[0]
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
        args = ""

        b, c, mess = BinCode()
        if mess != "":
            QMessageBox.critical(self,"Error", mess, QMessageBox.Ok, 0)
        else:
            args = c

            outfile = self.DialogCollector.ECSConversionDialog.resultFileName()

            args += " --no-write "
            args += " --case "
            args += os.path.join(thePath, outfile)
            args += " --post-volume "
            args += " med "
            args += CFDSTUDYGUI_DataModel._GetPath(sobj)

            log.debug("slotMeshConvertToMed -> args = %s" % args)
            dlg = CFDSTUDYGUI_CommandMgr.CFDSTUDYGUI_QProcessDialog(sgPyQt.getDesktop(),
                                                                    self.tr("STMSG_ECS_CONVERT"),
                                                                    [args],
                                                                    sobj.GetFather(),
                                                                    thePath)
            dlg.show()


    def slotCopyCaseFile(self):
        """
        Copy data xml file from a study case to another with popup menu attached to the xml file
        Copy into another case: COPY_CASE_FILE_ACTION_TEXT
        """
        sobj = self._singleSelectedObject()
        if sobj == None:
            return

        self.DialogCollector.CopyDialog.setCurrentObject(sobj)
        self.DialogCollector.CopyDialog.show()

        if not self.DialogCollector.CopyDialog.result() == QDialog.Accepted:
            return

        # update Object Browser
        # aDirPath: path directory where the xml file is copied
        aDirPath = self.DialogCollector.CopyDialog.destCaseName()
        aDirObject = CFDSTUDYGUI_DataModel.findMaxDeepObject(aDirPath)

        if aDirObject != None:
            self.updateObjBrowser(CFDSTUDYGUI_DataModel.GetCase(aDirObject))

        # BUG si je fais directement: self.updateObjBrowser(aDirObject)


    def slotCheckCompilation(self):
        """
        """
        #get current selection
        sobj = self._singleSelectedObject()
        if sobj == None:
            return

        # get current case
        aCase = CFDSTUDYGUI_DataModel.GetCase(sobj)
        if aCase == None:
            return

        cmd = self.__compile(aCase)
        if cmd != "":
            aChildList = CFDSTUDYGUI_DataModel.ScanChildren(aCase, "SRC")
            aSRCObj =  aChildList[0]
            aSRCPath = CFDSTUDYGUI_DataModel._GetPath(aSRCObj)

            dlg = CFDSTUDYGUI_CommandMgr.CFDSTUDYGUI_QProcessDialog(sgPyQt.getDesktop(),
                                                                    self.tr("STMSG_CHECK_COMPILATION"),
                                                                    [cmd],
                                                                    None,
                                                                    aSRCPath)
            dlg.show()


    def slotRunScript(self):
        """
        """
        log.debug("slotRunScript")
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

            dlg = CFDSTUDYGUI_CommandMgr.CFDSTUDYGUI_QProcessDialog(sgPyQt.getDesktop(),
                                                                    self.tr("STMSG_RUN_SCRIPT"),
                                                                    [path],
                                                                    father,
                                                                    fatherpath)
            dlg.show()


    def slotSaveDataFile(self):
        """
        Redirects B{Save} method to GUI of current solver
        """
        log.debug("slotSaveDataFile")
        if self._SolverGUI._CurrentWindow != None:
            if self._SolverGUI._CurrentWindow.case['xmlfile'] != "":
                self._SolverGUI._CurrentWindow.fileSave()
            else:
                self.slotSaveAsDataFile()


    def slotSaveAsDataFile(self):
        """
        Redirects B{SaveAs} method to GUI of current solver
        """
        log.debug("slotSaveAsDataFile")
        old_sobj = None
        new_sobj = None
        oldCase = self._SolverGUI.getCase(self._SolverGUI._CurrentWindow)
        oldStudy = CFDSTUDYGUI_DataModel.GetStudyByObj(oldCase)
        old_xml_file, xml_file = self._SolverGUI.SaveAsXmlFile()
        if old_xml_file == None and xml_file != None:
            #MP 25/04/2012 - A faire: tester si le fichier xml_file est deja ouvert dans une etude SALOME avec CFDSTUDYGUI_Management.py
            # classe CFDGUI_Management, methode findElem(xmlName, caseName, studyCFDName)
            # emettre un warning car on vient de sauvegarder dans un fichier xml existant et de plus ouvert dans une etude salome
            theNewStudyPath = os.path.dirname(os.path.dirname(os.path.dirname(xml_file)))
            study = CFDSTUDYGUI_DataModel.FindStudyByPath(theNewStudyPath)

            if study == None:
                theCaseName = os.path.basename(os.path.dirname(os.path.dirname(xml_file)))
                iok = CFDSTUDYGUI_DataModel._SetStudyLocation(theNewStudyPath, theCaseName)
                if iok:
                    study = CFDSTUDYGUI_DataModel.FindStudyByPath(theNewStudyPath)
                    obj = CFDSTUDYGUI_DataModel.checkPathUnderObject(study, xml_file)
                    if obj:
                        NewSObj = CFDSTUDYGUI_DataModel.getSObject(obj,os.path.basename(xml_file))
                        if  NewSObj != None:
                            self.OpenCFD_GUI(NewSObj)
                            self._SolverGUI.removeDockWindow(oldStudy.GetName(), oldCase.GetName(), "unnamed")
            else:
                theCaseName = os.path.basename(os.path.dirname(os.path.dirname(xml_file)))
                theCaseObj = CFDSTUDYGUI_DataModel.getSObject(study,theCaseName)
                if theCaseObj != None:
                    obj = CFDSTUDYGUI_DataModel.getSObject(theCaseObj,"DATA")
                    if obj != None:
                        NewSObj = CFDSTUDYGUI_DataModel.getSObject(obj,os.path.basename(xml_file))
                        if  NewSObj != None:
                            self._SolverGUI.removeDockWindow(oldStudy.GetName(), oldCase.GetName(), "unnamed")
                            self.CloseCFD_GUI(NewSObj)
                            self.OpenCFD_GUI(NewSObj)
                        else:
                            CFDSTUDYGUI_DataModel._CreateItem(obj,os.path.basename(xml_file))
                            NewSObj = CFDSTUDYGUI_DataModel.getSObject(obj,os.path.basename(xml_file))
                            if  NewSObj != None:
                                self._SolverGUI.removeDockWindow(study.GetName(),theCaseName , "unnamed")
                                self.OpenCFD_GUI(NewSObj)
                    else:
                        mess = "DATA directory is not found into Object Browser for case " +  theCaseName + "and study = " + study.GetName()
                        QMessageBox.critical(None, "Error", mess, QMessageBox.Ok, 0)
            return

        if xml_file != None and xml_file != old_xml_file and old_xml_file != None:
            theOldStudyPath = os.path.dirname(os.path.dirname(os.path.dirname(old_xml_file)))
            theOldStudyName = os.path.basename(theOldStudyPath)
            theNewStudyPath = os.path.dirname(os.path.dirname(os.path.dirname(xml_file)))
            theNewStudyName = os.path.basename(theNewStudyPath)
            oldStudy = CFDSTUDYGUI_DataModel.FindStudyByPath(theOldStudyPath)
            Old_obj = CFDSTUDYGUI_DataModel.checkPathUnderObject(oldStudy, old_xml_file) #parent DATA path object for old_xml_file
            if Old_obj:
                OldSobj = CFDSTUDYGUI_DataModel.getSObject(Old_obj,os.path.basename(old_xml_file))
                if OldSobj != None:
                    self.CloseCFD_GUI(OldSobj)

            if theOldStudyName == theNewStudyName:
                study = oldStudy
                obj = CFDSTUDYGUI_DataModel.checkPathUnderObject(study, xml_file) #parent DATA path object for xml_file
                if obj:
                    if os.path.exists(xml_file):
                        CFDSTUDYGUI_DataModel._CreateItem(obj,os.path.basename(xml_file))
                        NewSObj = CFDSTUDYGUI_DataModel.getSObject(obj,os.path.basename(xml_file))
                        if  NewSObj != None:
                            self.OpenCFD_GUI(NewSObj)

            else:
                study = CFDSTUDYGUI_DataModel.FindStudyByPath(theNewStudyPath)
                if study == None:
                    theCaseName = os.path.basename(os.path.dirname(os.path.dirname(xml_file)))
                    iok = CFDSTUDYGUI_DataModel._SetStudyLocation(theNewStudyPath, theCaseName)
                    if iok:
                        study = CFDSTUDYGUI_DataModel.FindStudyByPath(theNewStudyPath)
                        obj = CFDSTUDYGUI_DataModel.checkPathUnderObject(study, xml_file)
                        if obj:
                            NewSObj = CFDSTUDYGUI_DataModel.getSObject(obj,os.path.basename(xml_file))
                            if  NewSObj != None:
                                self.OpenCFD_GUI(NewSObj)
                else:
                   obj = CFDSTUDYGUI_DataModel.checkPathUnderObject(study, xml_file)
                   NewSObj = CFDSTUDYGUI_DataModel.getSObject(obj,os.path.basename(xml_file))
                   if  NewSObj != None:
                       self.CloseCFD_GUI(NewSObj)
                       self.OpenCFD_GUI(NewSObj)
                   else:
                       CFDSTUDYGUI_DataModel._CreateItem(obj,os.path.basename(xml_file))
                       NewSObj = CFDSTUDYGUI_DataModel.getSObject(obj,os.path.basename(xml_file))
                       if  NewSObj != None:
                           self.OpenCFD_GUI(NewSObj)


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


    def slotHelpLicense(self):
        self._SolverGUI.onSaturneHelpLicense()


    def slotHelpUserGuide(self):
        self._SolverGUI.onSaturneHelpManual()


    def slotHelpTutorial(self):
        self._SolverGUI.onSaturneHelpTutorial()


    def slotHelpTheory(self):
        self._SolverGUI.onSaturneHelpKernel()


    def slotHelpRefcard(self):
        self._SolverGUI.onSaturneHelpRefcard()


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
        elif theId in self._HelpActionIdMap:
            action_id = self._HelpActionIdMap[theId]

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
        elif theId in self._HelpActionIdMap:
            action_id = self._HelpActionIdMap[theId]

        if action_id == None:
            raise ActionError, "Invalid action id"

        return action_id


    def dskAgent(self):
        """
        Returns the dekstop Agent.
        """
        return self._DskAgent


    def disconnectSolverGUI(self):
        """
        Hide all the dock windows of CFDSTUDY, when activating another Salome Component
        We can have one or several of them with the right click on the main menu bar of
        Salome
        """
        log.debug("disconnectSolverGUI")
        #self._SolverGUI.disconnectDockWindowsBrowser()
        self._SolverGUI.disconnectDockWindows()


    def connectSolverGUI(self):
        """
        Show all the dock windows of CFDSTUDY, when activating another Salome Component
        """
        log.debug("connectSolverGUI")
        self._SolverGUI.connectDockWindows()
        #self._SolverGUI.connectDockWindowsBrowser()
