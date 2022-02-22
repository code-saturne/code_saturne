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
- BoundaryConditionsScalarView
"""

#-------------------------------------------------------------------------------
# Standard modules
#-------------------------------------------------------------------------------

import string, logging

#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------

from code_saturne.gui.base.QtCore    import *
from code_saturne.gui.base.QtGui     import *
from code_saturne.gui.base.QtWidgets import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.gui.case.BoundaryConditionsScalar import Ui_BoundaryConditionsScalar

from code_saturne.model.Common import GuiParam
from code_saturne.gui.base.QtPage import DoubleValidator, ComboModel, from_qvariant

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

        self.comboBoxScalarChoice.activated[str].connect(self.slotScalarTypeChoice)
        self.scalarChoiceModel = ComboModel(self.comboBoxScalarChoice, 1, 1)
        self.scalarChoiceModel.addItem(self.tr("Value"), 'dirichlet')
        self.scalarChoiceModel.addItem(self.tr("Flux"), 'flux')

        self.lineEditScalar.textChanged[str].connect(self.__slotScalar)

        validatorScalar = DoubleValidator(self.lineEditScalar)

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

        # For walls, current field is set to -1, hence the need for the list
        if boundary.getNature() == 'wall':
            fieldList = SpeciesModel(self.case).getFieldIdList(include_none=True)
        else:
            fieldList = [str(self.__currentField)]

        self.__scalarsList = []
        for f_id in fieldList:
            for species in SpeciesModel(self.case).getScalarByFieldId(f_id):
                self.__scalarsList.append((f_id, species))

        if len(self.__scalarsList) > 0 :
            for nb in range(len(self.__Scalarmodel.getItems())):
                self.__Scalarmodel.delItem(0)

            for var in self.__scalarsList :
                name = SpeciesModel(self.case).getScalarLabelByName(var[1])
                self.__Scalarmodel.addItem(self.tr(name), var[1])

            self.__currentScalar = self.__scalarsList[0]

            _f0, _s0 = self.__currentScalar

            self.__Scalarmodel.setItem(str_model=_s0)

            choice0 = self.__boundary.getScalarChoice(_f0, _s0)
            self.scalarChoiceModel.setItem(str_model=choice0)

            val = self.__boundary.getScalarValue(_f0, _s0)
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
        for _s in self.__scalarsList:
            if _s[1] == self.__Scalarmodel.dicoV2M[str(text)]:
                self.__currentScalar = _s
                break

        _f_id, _sname = self.__currentScalar

        # Update condition type choice
        scalarModel = self.__boundary.getScalarChoice(_f_id, _sname)
        self.scalarChoiceModel.setItem(str_model=scalarModel)

        # Update value
        val = self.__boundary.getScalarValue(_f_id, _sname)
        self.lineEditScalar.setText(str(val))


    @pyqtSlot(str)
    def slotScalarTypeChoice(self, text):
        """
        INPUT Condition type choice for scalar
        """
        choice = self.scalarChoiceModel.dicoV2M[str(text)]

        _f_id, _sname = self.__currentScalar

        self.__boundary.setScalarChoice(_f_id, _sname, choice)

    @pyqtSlot(str)
    def __slotScalar(self, text):
        """
        INPUT non condensable value
        """
        if self.lineEditScalar.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)

            _f_id, _sname = self.__currentScalar
            self.__boundary.setScalarValue(_f_id, _sname, value)


#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
