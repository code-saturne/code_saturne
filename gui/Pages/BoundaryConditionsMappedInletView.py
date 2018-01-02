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
This module contains the following classes:
- BoundaryConditionsMappedInletView
"""

#-------------------------------------------------------------------------------
# Standard modules
#-------------------------------------------------------------------------------

import string, logging

#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------

from code_saturne.Base.QtCore    import *
from code_saturne.Base.QtGui     import *
from code_saturne.Base.QtWidgets import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.Pages.BoundaryConditionsMappedInletForm import Ui_BoundaryConditionsMappedInletForm

from code_saturne.Base.Toolbox import GuiParam
from code_saturne.Base.QtPage import DoubleValidator, ComboModel, from_qvariant
from code_saturne.Pages.LocalizationModel import LocalizationModel, Zone
from code_saturne.Pages.Boundary import Boundary
from code_saturne.Pages.CompressibleModel import CompressibleModel

from code_saturne.Pages.QMeiEditorView import QMeiEditorView

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("BoundaryConditionsMappedInletView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class BoundaryConditionsMappedInletView(QWidget, Ui_BoundaryConditionsMappedInletForm):
    """
    Boundary condition for velocity in inlet, without particular physics.
    """
    def __init__(self, parent):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_BoundaryConditionsMappedInletForm.__init__(self)
        self.setupUi(self)


    def setup(self, case):
        """
        Setup the widget
        """
        self.__case = case
        self.__boundary = None

        self.__case.undoStopGlobal()

        # Connections
        self.groupBoxMappedInlet.clicked[bool].connect(self.__slotMappedInlet)

        self.lineEditTranslationX.textChanged[str].connect(self.__slotTrX)
        self.lineEditTranslationY.textChanged[str].connect(self.__slotTrY)
        self.lineEditTranslationZ.textChanged[str].connect(self.__slotTrZ)

        # Validators
        validatorX = DoubleValidator(self.lineEditTranslationX)
        validatorY = DoubleValidator(self.lineEditTranslationY)
        validatorZ = DoubleValidator(self.lineEditTranslationZ)

        # Apply validators
        self.lineEditTranslationX.setValidator(validatorX)
        self.lineEditTranslationY.setValidator(validatorY)
        self.lineEditTranslationZ.setValidator(validatorZ)

        self.__case.undoStartGlobal()


    def showWidget(self, boundary):
        """
        Show the widget
        """
        self.__boundary = boundary

        checked = False
        hide = False

        mdl = CompressibleModel(self.__case)
        if mdl.getCompressibleModel() != "off":
            hide = True

        if not hide and self.__boundary.getMappedInletStatus() == "on":
            checked = True

        if checked:
            self.groupBoxMappedInlet.setChecked(True)
        else:
            self.groupBoxMappedInlet.setChecked(False)
        self.__slotMappedInlet(checked)

        if hide:
            self.hide()
        else:
            self.show()


    def hideWidget(self):
        """
        Hide all
        """
        self.hide()


    @pyqtSlot(bool)
    def __slotMappedInlet(self, checked):
        """
        Private slot.

        Activates mapped inlet boundary condition.

        @type checked: C{True} or C{False}
        @param checked: if C{True}, shows the QGroupBox mapped inlet parameters.
        """
        self.groupBoxMappedInlet.setFlat(not checked)

        if checked:
            self.__boundary.setMappedInletStatus("on")
            self.labelDirection.show()
            self.frameDirectionCoordinates.show()
            x = self.__boundary.getMappedInletTranslation('translation_x')
            y = self.__boundary.getMappedInletTranslation('translation_y')
            z = self.__boundary.getMappedInletTranslation('translation_z')
            self.lineEditTranslationX.setText(str(x))
            self.lineEditTranslationY.setText(str(y))
            self.lineEditTranslationZ.setText(str(z))
        else:
            self.__boundary.setMappedInletStatus("off")
            self.labelDirection.hide()
            self.frameDirectionCoordinates.hide()


    @pyqtSlot(str)
    def __slotTrX(self, text):
        """
        INPUT value into direction of mapping translation
        """
        if self.lineEditTranslationX.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.__boundary.setMappedInletTranslation('translation_x', value)


    @pyqtSlot(str)
    def __slotTrY(self, text):
        """
        INPUT value into direction of mapping translation
        """
        if self.lineEditTranslationY.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.__boundary.setMappedInletTranslation('translation_y', value)


    @pyqtSlot(str)
    def __slotTrZ(self, text):
        """
        INPUT value into direction of mapping translation
        """
        if self.lineEditTranslationZ.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.__boundary.setMappedInletTranslation('translation_z', value)


    def tr(self, text):
        """
        Translation
        """
        return text

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
