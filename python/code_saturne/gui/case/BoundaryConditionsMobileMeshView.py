# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2022 EDF S.A.
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

from code_saturne.gui.base.QtCore    import *
from code_saturne.gui.base.QtGui     import *
from code_saturne.gui.base.QtWidgets import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.model.Common import GuiParam
from code_saturne.gui.base.QtPage import ComboModel

from code_saturne.gui.case.BoundaryConditionsMobileMeshForm import Ui_BoundaryConditionsMobileMeshForm
from code_saturne.model.MobileMeshModel import MobileMeshModel
from code_saturne.model.LocalizationModel import LocalizationModel, Zone
from code_saturne.model.Boundary import Boundary

from code_saturne.gui.case.QMegEditorView import QMegEditorView
from code_saturne.model.NotebookModel import NotebookModel

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
        self.case = case
        self.__boundary = None

        self.case.undoStopGlobal()

        self.__model = MobileMeshModel(self.case)
        self.notebook = NotebookModel(self.case)

        self.__comboModel = ComboModel(self.comboMobilBoundary, 6, 1)
        self.__comboModel.addItem(self.tr("Fixed boundary"), "fixed_boundary")
        self.__comboModel.addItem(self.tr("Sliding boundary"), "sliding_boundary")
        self.__comboModel.addItem(self.tr("Internal coupling"), "internal_coupling")
        self.__comboModel.addItem(self.tr("External coupling"), "external_coupling")
        self.__comboModel.addItem(self.tr("Fixed velocity"), "fixed_velocity")
        self.__comboModel.addItem(self.tr("Fixed displacement"), "fixed_displacement")

        self.comboMobilBoundary.activated[str].connect(self.__slotCombo)
        self.pushButtonMobilBoundary.clicked.connect(self.__slotFormula)

        self.case.undoStartGlobal()


    @pyqtSlot()
    def __slotFormula(self):
        """
        Run formula editor.
        """
        exp = self.__boundary.getALEFormula()
        c = self.__boundary.getALEChoice();

        if c == "fixed_velocity":
            if not exp:
                exp = 'mesh_velocity[0] = 0.;\nmesh_velocity[1] = 0.;\nmesh_velocity[2] = 0.;'
            req = [('mesh_velocity[0]', 'Fixed velocity of the mesh'),
                   ('mesh_velocity[1]', 'Fixed velocity of the mesh'),
                   ('mesh_velocity[2]', 'Fixed velocity of the mesh')]
            exa = 'mesh_velocity[0] = 0.;\nmesh_velocity[1] = 0.;\nmesh_velocity[2] = 1.;'
        elif c == "fixed_displacement":
            if not exp:
                exp = 'mesh_displacement[0] = 0.;\nmesh_displacement[1] = 0.;\nmesh_displacement[2] = 0.;'
            req = [('mesh_displacement[0]', 'Fixed displacement of the mesh'),
                   ('mesh_displacement[1]', 'Fixed displacement of the mesh'),
                   ('mesh_displacement[2]', 'Fixed displacement of the mesh')]
            exa = 'mesh_displacement[0] = 0.;\nmesh_displacement[1] = 0.;\nmesh_displacement[2] = 1.;'

        sym = [('x', "X face's gravity center"),
               ('y', "Y face's gravity center"),
               ('z', "Z face's gravity center"),
               ('dt', 'time step'),
               ('t', 'current time'),
               ('iter', 'number of iteration'),
               ('surface', 'Boundary zone surface')]

        for (nme, val) in self.notebook.getNotebookList():
            sym.append((nme, 'value (notebook) = ' + str(val)))

        dialog = QMegEditorView(parent = self,
                                function_type = 'bnd',
                                zone_name     = self.__boundary._label,
                                variable_name = 'mesh_velocity',
                                expression    = exp,
                                required      = req,
                                symbols       = sym,
                                condition     = c,
                                examples      = exa)

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
        exp = self.__boundary.getALEFormula()

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


#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
