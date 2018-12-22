# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2019 EDF S.A.
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
This module contains the following classes:
- BoundaryConditionsMobileMeshView
"""

#-------------------------------------------------------------------------------
# Standard modules
#-------------------------------------------------------------------------------

import logging

#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------

from code_saturne.Base.QtCore    import *
from code_saturne.Base.QtGui     import *
from code_saturne.Base.QtWidgets import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.Base.Toolbox import GuiParam
from code_saturne.Base.QtPage import ComboModel

from code_saturne.Pages.BoundaryConditionsMobileMeshForm import Ui_BoundaryConditionsMobileMeshForm
from code_saturne.Pages.MobileMeshModel import MobileMeshModel
from code_saturne.Pages.LocalizationModel import LocalizationModel, Zone
from code_saturne.Pages.Boundary import Boundary

from code_saturne.Pages.QMeiEditorView import QMeiEditorView
from code_saturne.Pages.NotebookModel import NotebookModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("BoundaryConditionsMobileMeshView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class BoundaryConditionsMobileMeshView(QWidget, Ui_BoundaryConditionsMobileMeshForm):
    """
    Boundary condifition for mobil mesh (ALE and/or Fluid-interaction)
    """
    def __init__(self, parent):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_BoundaryConditionsMobileMeshForm.__init__(self)
        self.setupUi(self)


    def setup(self, case):
        """
        Setup the widget
        """
        self.__case = case
        self.__boundary = None

        self.__case.undoStopGlobal()

        self.__model = MobileMeshModel(self.__case)
        self.notebook = NotebookModel(self.__case)

        self.__comboModel = ComboModel(self.comboMobilBoundary, 6, 1)
        self.__comboModel.addItem(self.tr("Fixed boundary"), "fixed_boundary")
        self.__comboModel.addItem(self.tr("Sliding boundary"), "sliding_boundary")
        self.__comboModel.addItem(self.tr("Internal coupling"), "internal_coupling")
        self.__comboModel.addItem(self.tr("External coupling"), "external_coupling")
        self.__comboModel.addItem(self.tr("Fixed velocity"), "fixed_velocity")
        self.__comboModel.addItem(self.tr("Fixed displacement"), "fixed_displacement")

        self.comboMobilBoundary.activated[str].connect(self.__slotCombo)
        self.pushButtonMobilBoundary.clicked.connect(self.__slotFormula)

        self.__case.undoStartGlobal()


    @pyqtSlot()
    def __slotFormula(self):
        """
        Run formula editor.
        """
        exp = self.__boundary.getFormula()
        aleChoice = self.__boundary.getALEChoice();

        if aleChoice == "fixed_velocity":
            if not exp:
                exp = 'mesh_velocity_U ='
            req = [('mesh_velocity_U', 'Fixed velocity of the mesh'),
                   ('mesh_velocity_V', 'Fixed velocity of the mesh'),
                   ('mesh_velocity_W', 'Fixed velocity of the mesh')]
            exa = 'mesh_velocity_U = 1000;\nmesh_velocity_V = 1000;\nmesh_velocity_W = 1000;'
        elif aleChoice == "fixed_displacement":
            if not exp:
                exp = 'mesh_x ='
            req = [('mesh_x', 'Fixed displacement of the mesh'),
                   ('mesh_y', 'Fixed displacement of the mesh'),
                   ('mesh_z', 'Fixed displacement of the mesh')]
            exa = 'mesh_x = 1000;\nmesh_y = 1000;\nmesh_z = 1000;'

        symbs = [('dt', 'time step'),
                 ('t', 'current time'),
                 ('iter', 'number of iteration')]

        for (nme, val) in self.notebook.getNotebookList():
            symbs.append((nme, 'value (notebook) = ' + str(val)))

        dialog = QMeiEditorView(self,
                                check_syntax = self.__case['package'].get_check_syntax(),
                                expression = exp,
                                required   = req,
                                symbols    = symbs,
                                examples   = exa)
        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormulaMobileMeshBoundary -> %s" % str(result))
            self.__boundary.setFormula(str(result))
            self.pushButtonMobilBoundary.setStyleSheet("background-color: green")
            self.pushButtonMobilBoundary.setToolTip(result)


    @pyqtSlot(str)
    def __slotCombo(self, text):
        """
        Called when the combobox changed.
        """
        modelData = self.__comboModel.dicoV2M[str(text)]

        if modelData == self.__boundary.getALEChoice():
            return

        self.__boundary.setALEChoice(modelData)
        exp = self.__boundary.getFormula()

        # Hide/Show formula button.
        # Formula is always reset when changing values, so set
        # color to red.
        if modelData in ["fixed_velocity", "fixed_displacement"]:
            self.pushButtonMobilBoundary.show()
        else:
            self.pushButtonMobilBoundary.hide()
        if exp:
            self.pushButtonMobilBoundary.setStyleSheet("background-color: red")
            self.pushButtonMobilBoundary.setToolTip(exp)
        else:
            self.pushButtonMobilBoundary.setStyleSheet("background-color: red")


    def showWidget(self, b):
        """
        Show the widget
        """
        if self.__model.getMethod() != "off":
            self.__boundary = b
            modelData = b.getALEChoice()
            self.__comboModel.setItem(str_model=modelData)
            if modelData in ["fixed_velocity", "fixed_displacement"]:
                self.pushButtonMobilBoundary.show()
            else:
                self.pushButtonMobilBoundary.hide()
            self.show()
        else:
            self.hideWidget()


    def hideWidget(self):
        """
        Hide all
        """
        self.hide()


    def tr(self, text):
        """
        Translation
        """
        return text

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
