# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2024 EDF S.A.
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
This module defines the 'Main fields initialization' page.

This module contains the following classes:
- ImmersedBoundariesFSIView
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import os, sys, string, types
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
from code_saturne.gui.base.QtPage import ComboModel, DoubleValidator
from code_saturne.gui.case.ImmersedBoundariesFSI import Ui_ImmersedBoundariesFSI
from code_saturne.model.ImmersedBoundariesModel import ImmersedBoundariesModel
from code_saturne.gui.case.QMegEditorView import QMegEditorView

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("ImmersedBoundariesFSIView")
log.setLevel(GuiParam.DEBUG)


#-------------------------------------------------------------------------------
# ImmersedBoundariesFSI class
#-------------------------------------------------------------------------------

class ImmersedBoundariesFSIView(QWidget, Ui_ImmersedBoundariesFSI):

    def __init__(self, parent=None):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_ImmersedBoundariesFSI.__init__(self)
        self.setupUi(self)

        self.case = None
        self.zone = None
        self.zone_id = None


    def setup(self, case, ibm, current_obj):
        self.case = case
        self.ibm = ibm
        self.current_obj = current_obj
        self.case.undoStopGlobal()

        #FSI velocity and position
        self.lineEditXInit.textChanged[str].connect(self.slotObjXinit)
        self.lineEditYInit.textChanged[str].connect(self.slotObjYinit)
        self.lineEditZInit.textChanged[str].connect(self.slotObjZinit)

        self.lineEditXEq.textChanged[str].connect(self.slotObjXeq)
        self.lineEditYEq.textChanged[str].connect(self.slotObjYeq)
        self.lineEditZEq.textChanged[str].connect(self.slotObjZeq)

        self.lineEditVelXInit.textChanged[str].connect(self.slotObjVelXinit)
        self.lineEditVelYInit.textChanged[str].connect(self.slotObjVelYinit)
        self.lineEditVelZInit.textChanged[str].connect(self.slotObjVelZinit)

        self.lineEditAccXInit.textChanged[str].connect(self.slotObjAccXinit)
        self.lineEditAccYInit.textChanged[str].connect(self.slotObjAccYinit)
        self.lineEditAccZInit.textChanged[str].connect(self.slotObjAccZinit)

        self.lineEditOmegaX.textChanged[str].connect(self.slotObjOmegaX)
        self.lineEditOmegaY.textChanged[str].connect(self.slotObjOmegaY)
        self.lineEditOmegaZ.textChanged[str].connect(self.slotObjOmegaZ)

        self.lineEditOmegaX.textChanged[str].connect(self.slotObjAnglesXinit)
        self.lineEditOmegaY.textChanged[str].connect(self.slotObjAnglesYinit)
        self.lineEditOmegaZ.textChanged[str].connect(self.slotObjAnglesZinit)

        # Connection to block displacement and rotation for FSI object
        self.checkBoxBlockDX.stateChanged.connect(self.slotCheckBoxBlockDX)
        self.checkBoxBlockDY.stateChanged.connect(self.slotCheckBoxBlockDY)
        self.checkBoxBlockDZ.stateChanged.connect(self.slotCheckBoxBlockDZ)
        self.checkBoxBlockRX.stateChanged.connect(self.slotCheckBoxBlockRX)
        self.checkBoxBlockRY.stateChanged.connect(self.slotCheckBoxBlockRY)
        self.checkBoxBlockRZ.stateChanged.connect(self.slotCheckBoxBlockRZ)

        # FSI parameters
        validatorMinCycle = DoubleValidator(self.lineEditMinSubCycles, min = 0.0)
        validatorMaxCycle = DoubleValidator(self.lineEditMaxSubCycles, min = 0.0)
        validatorCvCriteria = DoubleValidator(self.lineEditConvergenceCriteria, min = 0.0)

        validatorMinCycle.setExclusiveMin(True)
        validatorMaxCycle.setExclusiveMin(True)
        validatorCvCriteria.setExclusiveMin(True)

        self.lineEditMinSubCycles.setValidator(validatorMinCycle)
        self.lineEditMaxSubCycles.setValidator(validatorMaxCycle)
        self.lineEditConvergenceCriteria.setValidator(validatorCvCriteria)

        self.lineEditMinSubCycles.textChanged[str].connect(self.slotObjMinCycle)
        self.lineEditMaxSubCycles.textChanged[str].connect(self.slotObjMaxCycle)
        self.lineEditConvergenceCriteria.textChanged[str].connect(self.slotObjCvCriteria)

        self.update()
        self.case.undoStartGlobal()


    def update(self):

        if (self.ibm.getOnOff() == 'off' or self.ibm.getNumberOfObjects() == 0):
            return

        if self.ibm.getObjectFSI(self.current_obj) == 'on':
            self.groupBoxObjFSIProperties.show()
            self.groupBoxBlocking.show()
            self.groupBoxFSIParam.show()

            #FSI parameters
            self.lineEditMinSubCycles.setText(str(self.ibm.getObjectMinCycle(self.current_obj)))
            self.lineEditMaxSubCycles.setText(str(self.ibm.getObjectMaxCycle(self.current_obj)))
            self.lineEditConvergenceCriteria.setText(str(self.ibm.getObjectCvCriteria(self.current_obj)))

            #Blocking
            if self.ibm.getObjectBlockDX(self.current_obj) == 'off':
                self.checkBoxBlockDX.setChecked(False)
            else:
                self.checkBoxBlockDX.setChecked(True)

            if self.ibm.getObjectBlockDY(self.current_obj) == 'off':
                self.checkBoxBlockDY.setChecked(False)
            else:
                self.checkBoxBlockDY.setChecked(True)

            if self.ibm.getObjectBlockDZ(self.current_obj) == 'off':
                self.checkBoxBlockDZ.setChecked(False)
            else:
                self.checkBoxBlockDZ.setChecked(True)

            if self.ibm.getObjectBlockRX(self.current_obj) == 'off':
                self.checkBoxBlockRX.setChecked(False)
            else:
                self.checkBoxBlockRX.setChecked(True)

            if self.ibm.getObjectBlockRY(self.current_obj) == 'off':
                self.checkBoxBlockRY.setChecked(False)
            else:
                self.checkBoxBlockRY.setChecked(True)

            if self.ibm.getObjectBlockRZ(self.current_obj) == 'off':
                self.checkBoxBlockRZ.setChecked(False)
            else:
                self.checkBoxBlockRZ.setChecked(True)

            x0,y0,z0 = self.ibm.getObjectInitPosition(self.current_obj)
            self.lineEditXInit.setText(x0)
            self.lineEditYInit.setText(y0)
            self.lineEditZInit.setText(z0)

            xe,ye,ze = self.ibm.getObjectEqPosition(self.current_obj)
            self.lineEditXEq.setText(xe)
            self.lineEditYEq.setText(ye)
            self.lineEditZEq.setText(ze)

            vx,vy,vz = self.ibm.getObjectInitVel(self.current_obj)
            self.lineEditVelXInit.setText(vx)
            self.lineEditVelYInit.setText(vy)
            self.lineEditVelZInit.setText(vz)

            ax,ay,az = self.ibm.getObjectInitAcc(self.current_obj)
            self.lineEditAccXInit.setText(ax)
            self.lineEditAccYInit.setText(ay)
            self.lineEditAccZInit.setText(az)

            wx,wy,wz = self.ibm.getObjectAngularVel(self.current_obj)
            self.lineEditOmegaX.setText(wx)
            self.lineEditOmegaY.setText(wy)
            self.lineEditOmegaZ.setText(wz)

            thetax,thetay,thetaz = self.ibm.getObjectInitialAngles(self.current_obj)
            self.lineEditAnglesinitX.setText(thetax)
            self.lineEditAnglesinitY.setText(thetay)
            self.lineEditAnglesinitZ.setText(thetaz)

        else:
            self.comboBoxDensity.setEnabled(False)
            self.lineEditDensity.setReadOnly(True)
            self.lineEditDensity.setEnabled(False)

            self.comboBoxStiffness.setEnabled(False)
            self.lineEditStiffness.setReadOnly(True)
            self.lineEditStiffness.setEnabled(False)

            self.comboBoxDamping.setEnabled(False)
            self.lineEditDamping.setReadOnly(True)
            self.lineEditDamping.setEnabled(False)

            self.groupBoxObjFSIProperties.hide()
            self.groupBoxBlocking.hide()
            self.groupBoxFSIParam.hide()


    @pyqtSlot(str)
    def slotObjMinCycle(self, text):
        if self.lineEditMinSubCycles.validator().state == QValidator.Acceptable:
            val = int(text)
            self.ibm.setObjectMinCycle(self.current_obj, text)


    @pyqtSlot(str)
    def slotObjMaxCycle(self, text):
        if self.lineEditMaxSubCycles.validator().state == QValidator.Acceptable:
            val = int(text)
            self.ibm.setObjectMaxCycle(self.current_obj, text)


    @pyqtSlot(str)
    def slotObjCvCriteria(self, text):
        if self.lineEditConvergenceCriteria.validator().state == QValidator.Acceptable:
            val = float(text)
            self.ibm.setObjectCvCriteria(self.current_obj, text)


    @pyqtSlot(str)
    def slotObjXinit(self, text):
        if text == '':
            text=0.

        self.ibm.setObjectInitPosition(self.current_obj, xini=text)


    @pyqtSlot(str)
    def slotObjYinit(self, text):
        if text == '':
            text=0.

        self.ibm.setObjectInitPosition(self.current_obj, yini=text)


    @pyqtSlot(str)
    def slotObjZinit(self, text):
        if text == '':
            text=0.

        self.ibm.setObjectInitPosition(self.current_obj, zini=text)


    @pyqtSlot(str)
    def slotObjXeq(self, text):
        if text == '':
            text=0.

        self.ibm.setObjectEqPosition(self.current_obj, xeq=text)


    @pyqtSlot(str)
    def slotObjYeq(self, text):
        if text == '':
            text=0.

        self.ibm.setObjectEqPosition(self.current_obj, yeq=text)


    @pyqtSlot(str)
    def slotObjZeq(self, text):
        if text == '':
            text=0.

        self.ibm.setObjectEqPosition(self.current_obj, zeq=text)


    @pyqtSlot(str)
    def slotObjVelXinit(self, text):
        if text == '':
            text=0.

        self.ibm.setObjectInitVel(self.current_obj, vx=text)


    @pyqtSlot(str)
    def slotObjVelYinit(self, text):
        if text == '':
            text=0.

        self.ibm.setObjectInitVel(self.current_obj, vy=text)


    @pyqtSlot(str)
    def slotObjVelZinit(self, text):
        if text == '':
            text=0.

        self.ibm.setObjectInitVel(self.current_obj, vz=text)


    @pyqtSlot(str)
    def slotObjAccXinit(self, text):
        if text == '':
            text=0.

        self.ibm.setObjectInitAcc(self.current_obj, ax=text)

    @pyqtSlot(str)
    def slotObjAccYinit(self, text):
        if text == '':
            text=0.

        self.ibm.setObjectInitAcc(self.current_obj, ay=text)

    @pyqtSlot(str)
    def slotObjAccZinit(self, text):
        if text == '':
            text=0.

        self.ibm.setObjectInitAcc(self.current_obj, az=text)

    @pyqtSlot(str)
    def slotObjOmegaX(self, text):
        if text == '':
            text=0.

        self.ibm.setObjectAngularVel(self.current_obj, wx=text)

    @pyqtSlot(str)
    def slotObjOmegaY(self, text):
        if text == '':
            text=0.

        self.ibm.setObjectAngularVel(self.current_obj, wy=text)

    @pyqtSlot(str)
    def slotObjOmegaZ(self, text):
        if text == '':
            text=0

        self.ibm.setObjectAngularVel(self.current_obj, wz=text)

    @pyqtSlot(str)
    def slotObjAnglesXinit(self, text):
        if text == '':
            text=0.

        self.ibm.setObjectInitialAngles(self.current_obj, theta_x=text)

    @pyqtSlot(str)
    def slotObjAnglesYinit(self, text):
        if text == '':
            text=0.

        self.ibm.setObjectInitialAngles(self.current_obj, theta_y=text)

    @pyqtSlot(str)
    def slotObjAnglesZinit(self, text):
        if text == '':
            text=0.

        self.ibm.setObjectInitialAngles(self.current_obj, theta_z=text)

    @pyqtSlot()
    def slotFormulaRho(self): # TODO for FSI
        """
        User formula for density of the FSI object
        """
        objId = self.current_obj

        exp, req, sym  = self.ibm.getFSIFormulaRho(objId-1)
        exa = """rho = 1.8;"""

        name = self.ibm.getObjectName(objId)

        dialog = QMegEditorView(parent        = self,
                                function_type = 'ibm_vol',
                                zone_name     = name,
                                variable_name = 'porous_FSI_density',
                                expression    = exp,
                                required      = req,
                                symbols       = sym,
                                known_fields  = [],
                                examples      = exa)

        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormulaPorousRho -> %s" % str(result))
            self.pushButtonDensity.setStyleSheet("background-color: green")
            self.pushButtonDensity.setToolTip(exp)
            self.ibm.setFSIObjectRhoFormula(objId-1, result)


    @pyqtSlot(int)
    def slotCheckBoxBlockDX(self, val):
        text = ''
        if val == 0:
            text = 'off'
        else:
            text = 'on'

        self.ibm.setObjectBlockDX(self.current_obj, text)


    @pyqtSlot(int)
    def slotCheckBoxBlockDY(self, val):
        text = ''
        if val == 0:
            text = 'off'
        else:
            text = 'on'

        self.ibm.setObjectBlockDY(self.current_obj, text)


    @pyqtSlot(int)
    def slotCheckBoxBlockDZ(self, val):
        text = ''
        if val == 0:
            text = 'off'
        else:
            text = 'on'

        self.ibm.setObjectBlockDZ(self.current_obj, text)


    @pyqtSlot(int)
    def slotCheckBoxBlockRX(self, val):
        text = ''
        if val == 0:
            text = 'off'
        else:
            text = 'on'

        self.ibm.setObjectBlockRX(self.current_obj, text)


    @pyqtSlot(int)
    def slotCheckBoxBlockRY(self, val):
        text = ''
        if val == 0:
            text = 'off'
        else:
            text = 'on'

        self.ibm.setObjectBlockRY(self.current_obj, text)


    @pyqtSlot(int)
    def slotCheckBoxBlockRZ(self, val):
        text = ''
        if val == 0:
            text = 'off'
        else:
            text = 'on'

        self.ibm.setObjectBlockRZ(self.current_obj, text)
