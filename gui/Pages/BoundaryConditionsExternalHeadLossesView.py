# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2020 EDF S.A.
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
- BoundaryConditionsExternalHeadLossesView
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

from code_saturne.model.Common import GuiParam
from code_saturne.Base.QtPage import from_qvariant, to_text_string

from code_saturne.Pages.BoundaryConditionsExternalHeadLossesForm import Ui_BoundaryConditionsExternalHeadLossesForm

from code_saturne.model.LocalizationModel import LocalizationModel, Zone
from code_saturne.model.Boundary import Boundary
from code_saturne.Pages.QMegEditorView import QMegEditorView
from code_saturne.model.NotebookModel import NotebookModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("BoundaryConditionsExternalHeadLossesView")
log.setLevel(GuiParam.DEBUG)


#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class BoundaryConditionsExternalHeadLossesView(QWidget, Ui_BoundaryConditionsExternalHeadLossesForm):
    def __init__(self, parent):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_BoundaryConditionsExternalHeadLossesForm.__init__(self)
        self.setupUi(self)


    def setup(self, case):
        """
        Setup the widget
        """
        self.case = case
        self.__boundary = None
        self.notebook = NotebookModel(self.case)

        self.case.undoStopGlobal()

        self.pushButtonHeadLossesFormula.clicked.connect(self.slotHeadLossesFormula)

        self.case.undoStartGlobal()


    def showWidget(self, b):
        """
        Show the widget
        """
        label = b.getLabel()
        self.__boundary = Boundary('free_inlet_outlet', label, self.case)
        exp = self.__boundary.getHeadLossesFormula()
        if exp:
            self.pushButtonHeadLossesFormula.setStyleSheet("background-color: green")
            self.pushButtonHeadLossesFormula.setToolTip(exp)
        else:
            self.pushButtonHeadLossesFormula.setStyleSheet("background-color: red")

        self.show()


    def hideWidget(self):
        """
        Hide all
        """
        self.hide()


    @pyqtSlot()
    def slotHeadLossesFormula(self):
        """
        """
        exp = self.__boundary.getHeadLossesFormula()
        if not exp:
           exp = "K = 0.;"

        req = [('K', 'External head losses')]

        exa = "K = 0.;"

        sym = [('x', "X face's gravity center"),
               ('y', "Y face's gravity center"),
               ('z', "Z face's gravity center"),
               ('dt', 'time step'),
               ('t', 'current time'),
               ('iter', 'number of iteration'),
               ('surface', 'Boundary zone surface')]

        for (nme, val) in self.notebook.getNotebookList():
            sym.append((nme, 'value (notebook) = ' + str(val)))

        dialog = QMegEditorView(parent        = self,
                                function_type = 'bnd',
                                zone_name     = self.__boundary._label,
                                variable_name = 'head_loss',
                                expression    = exp,
                                required      = req,
                                symbols       = sym,
                                condition     = 'formula',
                                examples      = exa)

        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormulaDirection -> %s" % str(result))
            self.__boundary.setHeadLossesFormula(str(result))
            self.pushButtonHeadLossesFormula.setStyleSheet("background-color: green")
            self.pushButtonHeadLossesFormula.setToolTip(result)


    def tr(self, text):
        """
        Translation
        """
        return text

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
