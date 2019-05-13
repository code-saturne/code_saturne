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
- BoundaryConditionsNonCondensableView
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

from code_saturne.Pages.BoundaryConditionsNonCondensable import Ui_BoundaryConditionsNonCondensable

from code_saturne.model.Common import GuiParam
from code_saturne.Base.QtPage import DoubleValidator, ComboModel, from_qvariant

from code_saturne.model.NonCondensableModel import NonCondensableModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("BoundaryConditionsNonCondensableView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class BoundaryConditionsNonCondensableView(QWidget, Ui_BoundaryConditionsNonCondensable) :
    """
    Boundary condition for non condensable
    """
    def __init__(self, parent):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_BoundaryConditionsNonCondensable.__init__(self)
        self.setupUi(self)

        # Connections
        self.comboBoxNonCondensable.activated[str].connect(self.__slotChoiceNonCondensable)

        self.__NonCondensablemodel = ComboModel(self.comboBoxNonCondensable, 1, 1)

        self.lineEditNonCondensable.textChanged[str].connect(self.__slotNonCondensable)

        validatorNonCond = DoubleValidator(self.lineEditNonCondensable, min = 0.)
        validatorNonCond.setExclusiveMin(False)

        self.lineEditNonCondensable.setValidator(validatorNonCond)


    def setup(self, case, fieldId):
        """
        Setup the widget
        """
        self.case = case
        self.__boundary = None
        self.__currentField = fieldId


    def showWidget(self, boundary):
        """
        Show the widget
        """
        self.__boundary = boundary

        NonCondensableList = NonCondensableModel(self.case).getNonCondensableByFieldId(self.__currentField)

        if len(NonCondensableList) > 0 :
            for nb in range(len(self.__NonCondensablemodel.getItems())):
                self.__NonCondensablemodel.delItem(0)

            for var in NonCondensableList :
                name = NonCondensableModel(self.case).getNonCondLabel(var)
                self.__NonCondensablemodel.addItem(self.tr(name), var)

            self.__currentNonCondensable = NonCondensableList[0]
            self.__NonCondensablemodel.setItem(str_model=self.__currentNonCondensable)
            val = self.__boundary.getNonCondensableValue(self.__currentField, self.__currentNonCondensable)
            self.lineEditNonCondensable.setText(str(val))
            self.show()
        else :
            self.hideWidget()


    def hideWidget(self):
        """
        Hide the widget
        """
        self.hide()


    @pyqtSlot(str)
    def __slotChoiceNonCondensable(self, text):
        """
        INPUT choice of non condensable
        """
        self.__currentNonCondensable = self.__NonCondensablemodel.dicoV2M[str(text)]

        val = self.__boundary.getNonCondensableValue(self.__currentField, self.__currentNonCondensable)
        self.lineEditNonCondensable.setText(str(val))


    @pyqtSlot(str)
    def __slotNonCondensable(self, text):
        """
        INPUT non condensable value
        """
        if self.lineEditNonCondensable.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.__boundary.setNonCondensableValue(self.__currentField, self.__currentNonCondensable, value)


    def tr(self, text):
        """
        Translation
        """
        return text

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------

