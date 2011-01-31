# -*- coding: utf-8 -*-
#============================================================================
#
#     This file is part of CFDSTUDY the plug-in for Salome
#     of Code_Saturne CFD tool.
#
#     Copyright (C) 1998-2010 EDF S.A., France
#
#     contact: saturne-support@edf.fr
#
#     CFDSTUDY is free software; you can redistribute it
#     and/or modify it under the terms of the GNU General Public License
#     as published by the Free Software Foundation; either version 2 of
#     the License, or (at your option) any later version.
#
#     CFDSTUDY is distributed in the hope that it will be
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
#============================================================================

"""
CFDSTUDY
========

Main file of the CFD_STUDY module. Defines the standard SALOME callback.
"""

#-------------------------------------------------------------------------------
# Standard modules
#-------------------------------------------------------------------------------

import os
import logging

#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------

from PyQt4.QtGui import QApplication, QCursor, QDialog, QMessageBox, QDockWidget
from PyQt4.QtCore import Qt, QObject, QVariant, SIGNAL

#-------------------------------------------------------------------------------
# Salome modules
#-------------------------------------------------------------------------------

from SalomePyQt import WT_ObjectBrowser, WT_PyConsole, WT_LogWindow
import SALOMEDS
import SALOMEDS_Attributes_idl

#-------------------------------------------------------------------------------
# Application modules
#-------------------------------------------------------------------------------

import CFDSTUDYGUI_ActionsHandler
import CFDSTUDYGUI_DataModel
import CFDSTUDYGUI_Commons
import CFDSTUDYGUI_DesktopMgr

from CFDSTUDYGUI_Commons import CFD_Code, sg, sgPyQt
from CFDSTUDYGUI_Commons import CFD_Saturne, CFD_Neptune
from CFDSTUDYGUI_Commons import CheckCFD_CodeEnv

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("CFDSTUDYGUI")
log.setLevel(logging.DEBUG)
#log.setLevel(logging.NOTSET)

#-------------------------------------------------------------------------------
# Global definitions
#-------------------------------------------------------------------------------

studyId = 1
d_activation = {}

# Desktop manager: instance of CFDSTUDYGUI_DesktopMgr class to store the SALOME Workspace
_DesktopMgr = CFDSTUDYGUI_DesktopMgr.CFDSTUDYGUI_DesktopMgr()

# ObjectTR is a convenient object for traduction purpose
ObjectTR = QObject()

#-------------------------------------------------------------------------------
# Callback GUI functions
#-------------------------------------------------------------------------------

def initialize():
    """
    This method is called when GUI module is being created and initialized.
    """
    # nothing to do here.
    log.debug("initialize")
    pass


def windows():
    """
    This method is called when GUI module is being created
    and initialized.
    Should return a layout of the SALOME dockable windows id's
    needed to be opened when module is activated.

    @return: layout of the SALOME dockable windows
    @rtype: C{Dictionary}
    """
    log.debug("windows")

    winMap = {}
    winMap[ WT_ObjectBrowser ] = Qt.LeftDockWidgetArea
    winMap[ WT_PyConsole ]     = Qt.BottomDockWidgetArea
    winMap[ WT_LogWindow ]     = Qt.BottomDockWidgetArea

    return winMap


def views():
    """
    This method is called when GUI module is being created and initialized.
    Should return a list of the SALOME view window types
    needed to be opened when module is activated.
    cf. SALOME_PYQT_Module.cxx : PyObjWrapper

    @return: list of the SALOME view window types
    @rtype: C{String} or C{list} of C{String}
    """
    log.debug("views")
    winList = "VTKViewer"
    #winList = "OCCViewer"
    #winList = "Plot2d"

    #winList = ""
    return winList


def setWorkSpace(ws):
    """
    Stores the SALOME Workspace I{ws} into the Desktop Manager.

    @type ws: C{QWidget}
    @param ws: main window's central widget
    """
    log.debug("setWorkSpace")

    dsk = sgPyQt.getDesktop()
    _DesktopMgr.setWorkspace(dsk, ws)


def createPreferences():
    """
    Manages the preferences QDialog of the module.
    """
    log.debug("createPreferences")
    sgPyQt = CFDSTUDYGUI_DataModel.sgPyQt
    tabId = sgPyQt.addPreference(ObjectTR.tr("CFDSTUDY_PREF_TAB"))
    genGroup = sgPyQt.addPreference(ObjectTR.tr("CFDSTUDY_PREF_GEN_GROUP"), tabId)
    sgPyQt.setPreferenceProperty(genGroup, "columns", QVariant(1))

    sgPyQt.addPreference(ObjectTR.tr("CFDSTUDY_PREF_EDITOR"), genGroup, 3, "CFDSTUDY", "Editor")
    sgPyQt.addPreference(ObjectTR.tr("CFDSTUDY_PREF_READER"), genGroup, 3, "CFDSTUDY", "Reader")
    sgPyQt.addPreference(ObjectTR.tr("CFDSTUDY_PREF_DIFF"  ), genGroup, 3, "CFDSTUDY", "Tool_for_diff")
    sgPyQt.addPreference(ObjectTR.tr("CFDSTUDY_PREF_DISPLAY"), genGroup, 3, "CFDSTUDY", "Display")

    genGroup = sgPyQt.addPreference(ObjectTR.tr("CFDSTUDY_PREF_ENV_GROUP"), tabId)
    sgPyQt.setPreferenceProperty(genGroup, "columns", QVariant(1))


def activate():
    """
    This method is called when GUI module is being activated.

    @rtype: C{True} or C{False}
    @return: C{True} only if the activation is successful.
    """
    log.debug("activate")
    global d_activation, studyId

    dsk = sgPyQt.getDesktop()
    studyId = sgPyQt.getStudyId()
    # instance of the CFDSTUDYGUI_ActionsHandler class for the current desktop
    ActionHandler = _DesktopMgr.getActionHandler(dsk)

    if studyId not in d_activation.keys():
        d_activation[studyId] = 1

    if d_activation[studyId] == 1:
	d_activation[studyId] = 0
        env_saturne = CheckCFD_CodeEnv(CFD_Saturne)
        env_neptune = CheckCFD_CodeEnv(CFD_Neptune)

        #log.debug("activate -> env_saturne = %s" % env_saturne)
        #log.debug("activate -> env_neptune = %s" % env_neptune)

        if not env_saturne and not env_neptune:
            mess = ObjectTR.tr("CFDSTUDY_INVALID_ENV")
            QMessageBox.critical(ActionHandler.dskAgent().workspace(),
                                 "Error", mess, QMessageBox.Ok, 0)

	    d_activation[studyId] = 1
            return False
        elif env_saturne:
            ActionHandler.DialogCollector.InfoDialog.setCode(CFD_Saturne, True)
        elif env_neptune:
            ActionHandler.DialogCollector.InfoDialog.setCode(CFD_Neptune, True)

        ActionHandler.DialogCollector.InfoDialog.exec_()

        if not ActionHandler.DialogCollector.InfoDialog.result() == QDialog.Accepted:
	    d_activation[studyId] = 1
            return False

    ActionHandler.connect(ActionHandler._SalomeSelection,
                          SIGNAL('currentSelectionChanged()'),
                          ActionHandler.updateActions)

    ActionHandler.connectSolverGUI()
    ActionHandler.updateObjBrowser()

    # Hide the Python Console window layout
    dsk = sgPyQt.getDesktop()
    ldockWindows = dsk.findChildren(QDockWidget)
    for dock in ldockWindows:
        dockTitle = dock.windowTitle()
        if str(dockTitle) == "Python Console":
            dock.setVisible(False)

    return True


def setSettings():
    """
    Stores the selected CFD code and updates action according with current
    selection and study states in the dekstop manager.
    """
    log.debug("setSettings")

    dsk = sgPyQt.getDesktop()
    ActionHandler = _DesktopMgr.getActionHandler(dsk)
    ActionHandler.onCFDCode()
    ActionHandler.updateActions()


def deactivate():
    """
    This method is called when GUI module is being deactivated.
    """
    log.debug("deactivate")
    dsk = sgPyQt.getDesktop()
    ActionHandler = _DesktopMgr.getActionHandler(dsk)
    ActionHandler.disconnect(ActionHandler._SalomeSelection, SIGNAL('currentSelectionChanged()'), ActionHandler.updateActions)
    ActionHandler.disconnectSolverGUI()


def createPopupMenu(popup, context):
    """
    This method is called when popup menu is requested by the user (right click).
    Should analyze the selection and fill in the popup menu with the corresponding actions.

    @type popup: C{QPopupMenu}
    @param popup: popup menu from the Object Browser.
    @type context: C{String}
    @param context: equal to 'ObjectBrowser' or 'VTKViewer' for example.
    """
    log.debug("createPopupMenu -> context = %s" % context)

    study = CFDSTUDYGUI_DataModel._getStudy()
    dsk = sgPyQt.getDesktop()
    ActionHandler = _DesktopMgr.getActionHandler(dsk)

    log.debug("createPopupMenu -> SelectedCount = %s" % sg.SelectedCount())

    id_flag = 0

    if sg.SelectedCount() <= 0:
        return
    elif sg.SelectedCount() == 1:
        entry = sg.getSelected(0)
        if entry != '':
            sobj = study.FindObjectID(entry)
            if sobj is not None:
                test, anAttr = sobj.FindAttribute("AttributeLocalID")
                if test:
                    id = anAttr._narrow(SALOMEDS.AttributeLocalID).Value()
                    if id >= 0:
                        ActionHandler.customPopup(id, popup)
                        if CFDSTUDYGUI_DataModel.isLinkPathObject(sobj):
                            popup.removeAction(ActionHandler.commonAction(CFDSTUDYGUI_ActionsHandler.EditAction))
                            popup.removeAction(ActionHandler.commonAction(CFDSTUDYGUI_ActionsHandler.MoveToDRAFTAction))
                            popup.removeAction(ActionHandler.commonAction(CFDSTUDYGUI_ActionsHandler.CopyCaseFileAction))
                    id_flag = id

    else:
#        flag = True
#        index = 0

#        #check for Pre MED files multi selection
#        while index < sg.SelectedCount() and flag :
#            sobj = study.FindObjectID(sg.getSelected(index))
#            flag = CFDSTUDYGUI_DataModel.checkPreMEDType(sobj)
#            index += 1
#
#        if not flag:
#            return
#
#        if flag:
#            #add MED conversion to popup
#            popup.addAction(ActionHandler.commonAction(CFDSTUDYGUI_ActionsHandler.ECSConvertAction))
#            return

        #add Display Mesh group to popup

        for i in range(sg.SelectedCount()):
            entry = sg.getSelected(i)
            if entry != '':
                sobj = study.FindObjectID(entry)
                if sobj is not None:
                    test, anAttr = sobj.FindAttribute("AttributeLocalID")
                    if test:
                        id = anAttr._narrow(SALOMEDS.AttributeLocalID).Value()
                        if id != CFDSTUDYGUI_DataModel.dict_object["PublishedMeshGroupIntoObjectBrowser"]:
                            id_flag = 0
                            break
                        else:
                            id_flag = id

            if id_flag == CFDSTUDYGUI_DataModel.dict_object["PublishedMeshGroupIntoObjectBrowser"]:
                ActionHandler.customPopup(id, popup)
                popup.removeAction(ActionHandler.commonAction(CFDSTUDYGUI_ActionsHandler.DisplayOnlyGroupMESHAction))

    if context == "VTKViewer":
        if id_flag == CFDSTUDYGUI_DataModel.dict_object["PublishedMeshGroupIntoObjectBrowser"]:
            ActionHandler.customPopup("VTKViewer", popup)





# TODELETE
#def processMgr():
#    """
#    Stores the SALOME desktop (i.e. C{QMainWindow}) in the dekstop manager.
#    """
#    log.debug("processMgr")
#    dsk = sgPyQt.getDesktop()
#    return _DesktopMgr.getProcessMgr(dsk)
