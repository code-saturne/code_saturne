# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2015 EDF S.A.
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
This module defines the values of reference.

This module contains the following classes and function:
- DarcyView
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import logging

#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------

from PyQt4.QtCore import *
from PyQt4.QtGui  import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.Base.Toolbox     import GuiParam
from code_saturne.Base.QtPage      import ComboModel, from_qvariant, DoubleValidator
from code_saturne.Pages.DarcyForm  import Ui_DarcyForm
from code_saturne.Pages.DarcyModel import DarcyModel

from code_saturne.Pages.QMeiEditorView import QMeiEditorView

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("DarcyView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class DarcyView(QWidget, Ui_DarcyForm):
    """
    Class to open Page.
    """
    def __init__(self, parent, case):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_DarcyForm.__init__(self)
        self.setupUi(self)

        self.case = case
        self.case.undoStopGlobal()
        self.mdl = DarcyModel(self.case)

        # ComboBox
        self.modelPermeability = ComboModel(self.comboBoxPermeability,2,1)
        self.modelDiffusion = ComboModel(self.comboBoxDiffusion,2,1)
        self.modelFlowType = ComboModel(self.comboBoxFlowType,2,1)
        self.modelCriterion = ComboModel(self.comboBoxCriterion,2,1)

        self.lineEditGx.setValidator(DoubleValidator(self.lineEditGx))
        self.lineEditGy.setValidator(DoubleValidator(self.lineEditGy))
        self.lineEditGz.setValidator(DoubleValidator(self.lineEditGz))

        self.modelPermeability.addItem(self.tr("isotropic"), 'isotropic')
        self.modelPermeability.addItem(self.tr("anisotropic"), 'anisotropic')
        self.modelDiffusion.addItem(self.tr("isotropic"), 'isotropic')
        self.modelDiffusion.addItem(self.tr("anisotropic"), 'anisotropic')
        self.modelFlowType.addItem(self.tr("steady"), 'steady')
        self.modelFlowType.addItem(self.tr("unsteady"), 'unsteady')
        self.modelCriterion.addItem(self.tr("over pressure"), 'pressure')
        self.modelCriterion.addItem(self.tr("over velocity"), 'velocity')

        # Connections
        self.connect(self.comboBoxPermeability, SIGNAL("activated(const QString&)"), self.slotPermeabilityType)
        self.connect(self.comboBoxDiffusion, SIGNAL("activated(const QString&)"), self.slotDiffusionType)
        self.connect(self.comboBoxFlowType, SIGNAL("activated(const QString&)"), self.slotFlowType)
        self.connect(self.comboBoxCriterion, SIGNAL("activated(const QString&)"), self.slotCriterion)
        self.connect(self.checkBoxGravity, SIGNAL("clicked()"), self.slotGravity)
        self.connect(self.lineEditGx, SIGNAL("textChanged(const QString &)"), self.slotGravityX)
        self.connect(self.lineEditGy, SIGNAL("textChanged(const QString &)"), self.slotGravityY)
        self.connect(self.lineEditGz, SIGNAL("textChanged(const QString &)"), self.slotGravityZ)

        self.initializeWidget()

        self.case.undoStartGlobal()


    @pyqtSignature("")
    def initializeWidget(self):
        """
        """
        value = self.mdl.getPermeabilityType()
        self.modelPermeability.setItem(str_model=value)

        value = self.mdl.getDiffusionType()
        self.modelDiffusion.setItem(str_model=value)

        value = self.mdl.getFlowType()
        self.modelFlowType.setItem(str_model=value)

        value = self.mdl.getCriterion()
        self.modelCriterion.setItem(str_model=value)

        if self.mdl.getGravity() == 'on':
            self.checkBoxGravity.setChecked(True)
            self.groupBoxGravity.show()
            gx, gy, gz = self.mdl.getGravityVector()
            self.lineEditGx.setText(str(gx))
            self.lineEditGy.setText(str(gy))
            self.lineEditGz.setText(str(gz))
        else:
            self.checkBoxGravity.setChecked(False)
            self.groupBoxGravity.hide()


    @pyqtSignature("const QString&")
    def slotPermeabilityType(self, text):
        """
        Input permeability type : isotrop or anisotrop.
        """
        mdl = self.modelPermeability.dicoV2M[str(text)]
        self.mdl.setPermeabilityType(mdl)


    @pyqtSignature("const QString&")
    def slotDiffusionType(self, text):
        """
        Input viscosity type : isotrop or anisotrop.
        """
        mdl = self.modelDiffusion.dicoV2M[str(text)]
        self.mdl.setDiffusionType(mdl)


    @pyqtSignature("const QString&")
    def slotFlowType(self, text):
        """
        Input flow type : steady or unsteady.
        """
        mdl = self.modelFlowType.dicoV2M[str(text)]
        self.mdl.setFlowType(mdl)


    @pyqtSignature("const QString&")
    def slotCriterion(self, text):
        """
        Input convergence criterion of Newton scheme.
        """
        mdl = self.modelCriterion.dicoV2M[str(text)]
        self.mdl.setCriterion(mdl)


    @pyqtSignature("")
    def slotGravity(self):
        """
        Input if gravity is taken into account or not
        """
        if self.checkBoxGravity.isChecked():
            self.mdl.setGravity('on')
            self.groupBoxGravity.show()
            gx, gy, gz = self.mdl.getGravityVector()
            self.lineEditGx.setText(str(gx))
            self.lineEditGy.setText(str(gy))
            self.lineEditGz.setText(str(gz))
        else:
            self.mdl.setGravity('off')


    @pyqtSignature("const QString&")
    def slotGravityX(self, text):
        """
        Input component X for gravity vector
        """
        if self.sender().validator().state == QValidator.Acceptable:
            val = float(text)
            self.mdl.setGravityVector("axis_x", val)


    @pyqtSignature("const QString&")
    def slotGravityY(self, text):
        """
        Input component Y for gravity vector
        """
        if self.sender().validator().state == QValidator.Acceptable:
            val = float(text)
            self.mdl.setGravityVector("axis_y", val)


    @pyqtSignature("const QString&")
    def slotGravityZ(self, text):
        """
        Input component Z for gravity vector
        """
        if self.sender().validator().state == QValidator.Acceptable:
            val = float(text)
            self.mdl.setGravityVector("axis_z", val)


    def tr(self, text):
        """
        Translation
        """
        return text

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
