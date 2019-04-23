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
- BoundaryConditionsScalarView
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

from code_saturne.Pages.BoundaryConditionsScalar import Ui_BoundaryConditionsScalar

from code_saturne.model.Common import GuiParam
from code_saturne.Base.QtPage import DoubleValidator, ComboModel, from_qvariant

from code_saturne.model.SpeciesModel import SpeciesModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("BoundaryConditionsSCalarView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class BoundaryConditionsScalarView(QWidget, Ui_BoundaryConditionsScalar) :
    """
    Boundary condition for non condensable
    """
    def __init__(self, parent):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_BoundaryConditionsScalar.__init__(self)
        self.setupUi(self)

        # Connections
        self.comboBoxScalar.activated[str].connect(self.__slotChoiceScalar)

        self.__Scalarmodel = ComboModel(self.comboBoxScalar, 1, 1)

        self.lineEditScalar.textChanged[str].connect(self.__slotScalar)

        validatorScalar = DoubleValidator(self.lineEditScalar, min = 0.)
        validatorScalar.setExclusiveMin(False)

        self.lineEditScalar.setValidator(validatorScalar)


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

        ScalarList = SpeciesModel(self.case).getScalarByFieldId(self.__currentField)

        if len(ScalarList) > 0 :
            for nb in range(len(self.__Scalarmodel.getItems())):
                self.__Scalarmodel.delItem(0)

            for var in ScalarList :
                name = SpeciesModel(self.case).getScalarLabelByName(var)
                self.__Scalarmodel.addItem(self.tr(name), var)
            self.__currentScalar = ScalarList[0]
            self.__Scalarmodel.setItem(str_model=self.__currentScalar)
            val = self.__boundary.getScalarValue(self.__currentField, self.__currentScalar)
            self.lineEditScalar.setText(str(val))
            self.show()
        else :
            self.hideWidget()


    def hideWidget(self):
        """
        Hide the widget
        """
        self.hide()


    @pyqtSlot(str)
    def __slotChoiceScalar(self, text):
        """
        INPUT choice of non condensable
        """
        self.__currentScalar = self.__Scalarmodel.dicoV2M[str(text)]

        val = self.__boundary.getScalarValue(self.__currentField, self.__currentScalar)
        self.lineEditScalar.setText(str(val))


    @pyqtSlot(str)
    def __slotScalar(self, text):
        """
        INPUT non condensable value
        """
        if self.lineEditScalar.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.__boundary.setScalarValue(self.__currentField, self.__currentScalar, value)


    def tr(self, text):
        """
        Translation
        """
        return text

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
