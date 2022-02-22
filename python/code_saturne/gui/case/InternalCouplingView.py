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
This module defines the internal coupling view data management.
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

from code_saturne.model.InternalCouplingModel import InternalCouplingModel
from code_saturne.gui.base.QtPage import ComboModel

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class InternalCouplingView(QWidget):

    def __init__(self, parent=None):
        """
        Constructor
        """

        QWidget.__init__(self, parent)

        self.initUi()

    def initUi(self):

        # Main layout
        vbox = QVBoxLayout()

        # Line with Combobox + checkbox
        hbox = QHBoxLayout()

        self.comboBoxScalars = QComboBox()
        self.comboBoxScalars.setObjectName("comboBoxScalars")

        self.comboBoxScalars.currentIndexChanged[str].connect(self.slotScalarName)

        self.checkBoxScalars = QCheckBox()
        self.checkBoxScalars.setObjectName("checkBoxScalars")
        self.checkBoxScalars.stateChanged.connect(self.slotScalarState)
        self.checkBoxScalars.setText("couple scalar in solid zones")

        spacerItem = QSpacerItem(40, 20, QSizePolicy.Expanding, QSizePolicy.Minimum)

        hbox.addWidget(self.comboBoxScalars)
        hbox.addItem(QSpacerItem(40, 2, QSizePolicy.Minimum, QSizePolicy.Minimum))
        hbox.addWidget(self.checkBoxScalars)
        hbox.addItem(spacerItem)

        vbox.addItem(hbox)

        vbox.addItem(QSpacerItem(20, 0, QSizePolicy.Minimum, QSizePolicy.Expanding))

        self.setLayout(vbox)


    def setup(self, case):

        self.case = case
        self.mdl = InternalCouplingModel(case)

        self.sca = None

        self.combo = ComboModel(self.comboBoxScalars, 1, 1)
        scalars_list = self.mdl.getListOfAvailableScalars()
        if scalars_list:
            for sca in scalars_list:
                self.combo.addItem(sca)
            self.has_scalars = True
            self.sca = scalars_list[0]
        else:
            self.sca = ""
            self.combo.addItem("")
            self.has_scalars = False
            self.checkBoxScalars.setChecked(False)

        self.comboBoxScalars.setEnabled(self.has_scalars)
        self.checkBoxScalars.setEnabled(self.has_scalars)

        if self.has_scalars:
            if self.mdl.isScalarCoupled(self.sca):
                self.checkBoxScalars.setChecked(True)


    @pyqtSlot(int)
    def slotScalarState(self, val):

        if val == 0:
            self.mdl.removeScalar(self.sca)
        else:
            self.mdl.addScalar(self.sca)

    @pyqtSlot(str)
    def slotScalarName(self, val):
        self.sca = str(val)

        scalar_is_coupled = self.mdl.isScalarCoupled(self.sca)

        self.checkBoxScalars.setChecked(scalar_is_coupled)

