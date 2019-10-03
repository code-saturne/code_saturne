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
This module contains the following classes and function:
- LabelDelegate
- LagrangianAdvancedOptionsDialogForm
- StandardItemModelCoals
- LagrangianView
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
from code_saturne.Base.QtPage import ComboModel, IntValidator, DoubleValidator
from code_saturne.Base.QtPage import from_qvariant, to_text_string
from code_saturne.Pages.LagrangianForm import Ui_LagrangianForm
from code_saturne.Pages.LagrangianAdvancedOptionsDialogForm import Ui_LagrangianAdvancedOptionsDialogForm
from code_saturne.model.LagrangianModel import LagrangianModel
from code_saturne.model.CoalCombustionModel import CoalCombustionModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("LagrangianView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Popup advanced options
#-------------------------------------------------------------------------------

class LagrangianAdvancedOptionsDialogView(QDialog, Ui_LagrangianAdvancedOptionsDialogForm):
    """
    Advanced dialog
    """
    def __init__(self, parent, case, default):
        """
        Constructor
        """
        QDialog.__init__(self, parent)

        Ui_LagrangianAdvancedOptionsDialogForm.__init__(self)
        self.setupUi(self)

        self.case = case
        self.case.undoStopGlobal()

        self.setWindowTitle(self.tr("Advanced options"))
        self.default = default
        self.result  = self.default.copy()

        # Combo model
        self.modelNORDRE = ComboModel(self.comboBoxNORDRE, 2, 1)
        self.modelNORDRE.addItem(self.tr("first-order scheme"),  "1")
        self.modelNORDRE.addItem(self.tr("second-order scheme"), "2")

        self.modelIDIRLA = ComboModel(self.comboBoxIDIRLA, 3, 1)
        self.modelIDIRLA.addItem(self.tr("X"), "1")
        self.modelIDIRLA.addItem(self.tr("Y"), "2")
        self.modelIDIRLA.addItem(self.tr("Z"), "3")
        self.modelIDIRLA.addItem(self.tr("Local projection"), "4")

        # Connections
        self.comboBoxNORDRE.activated[str].connect(self.slotNORDRE)
        self.checkBoxIDISTU.clicked.connect(self.slotIDISTU)
        self.checkBoxIDIFFL.clicked.connect(self.slotIDIFFL)
        self.groupBoxModel.clicked[bool].connect(self.slotModel)
        self.lineEditMODCPL.textChanged[str].connect(self.slotMODCPL)
        self.comboBoxIDIRLA.activated[str].connect(self.slotIDIRLA)

        validatorMODCPL = IntValidator(self.lineEditMODCPL, min=1)
        self.lineEditMODCPL.setValidator(validatorMODCPL)

        # initialize Widgets
        order = str(self.result['scheme_order'])
        self.modelNORDRE.setItem(str_model=order)

        if self.result['turbulent_dispertion'] == "on":
            self.checkBoxIDISTU.setChecked(True)
        else:
            self.checkBoxIDISTU.setChecked(False)

        if self.result['fluid_particles_turbulent_diffusion'] == "on":
            self.checkBoxIDIFFL.setChecked(True)
        else:
            self.checkBoxIDIFFL.setChecked(False)

        value = self.result['complete_model_iteration']
        if value > 0:
            self.lineEditMODCPL.setText(str(value))

            direction = self.result['complete_model_direction']
            self.modelIDIRLA.setItem(str_model=str(direction))
        else:
            self.groupBoxModel.setChecked(False)

        self.case.undoStartGlobal()


    @pyqtSlot(str)
    def slotNORDRE(self, text):
        """
        Input NORDRE.
        """
        value = self.modelNORDRE.dicoV2M[str(text)]
        self.result['scheme_order'] = value


    @pyqtSlot()
    def slotIDISTU(self):
        """
        Input IDISTU.
        """
        if self.checkBoxIDISTU.isChecked():
            status = "on"
        else:
            status = "off"
        self.result['turbulent_dispertion'] = status


    @pyqtSlot()
    def slotIDIFFL(self):
        """
        Input IDIFFL.
        """
        if self.checkBoxIDIFFL.isChecked():
            status = "on"
        else:
            status = "off"
        self.result['fluid_particles_turbulent_diffusion'] = status


    @pyqtSlot(bool)
    def slotModel(self, checked):
        if checked:
             value = self.default['complete_model_iteration']
             if value == 0:
                 value = 1
             self.result['complete_model_iteration'] = value
             self.lineEditMODCPL.setText(str(value))
        else:
             self.result['complete_model_iteration'] = 0


    @pyqtSlot(str)
    def slotMODCPL(self, text):
        """
        Input MODCPL.
        """
        if self.lineEditMODCPL.validator().state == QValidator.Acceptable:
            self.result['complete_model_iteration'] = from_qvariant(text, int)


    @pyqtSlot(str)
    def slotIDIRLA(self, text):
        """
        Input IDIRLA.
        """
        value = self.modelIDIRLA.dicoV2M[str(text)]
        self.result['complete_model_direction'] = value


    def get_result(self):
        """
        Method to get the result
        """
        return self.result


    def accept(self):
        """
        Method called when user clicks 'OK'
        """
        QDialog.accept(self)


    def reject(self):
        """
        Method called when user clicks 'Cancel'
        """
        self.result = self.default.copy()
        QDialog.reject(self)


    def tr(self, text):
        """
        Translation
        """
        return text
#-------------------------------------------------------------------------------
# Line edit delegate for the label
#-------------------------------------------------------------------------------

class LabelDelegate(QItemDelegate):
    """
    Use of a QLineEdit in the table.
    """
    def __init__(self, parent=None):
        QItemDelegate.__init__(self, parent)
        self.parent = parent
        self.old_plabel = ""


    def createEditor(self, parent, option, index):
        editor = QLineEdit(parent)
        self.old_label = ""
        rx = "[_a-zA-Z][_A-Za-z0-9]{1," + str(LABEL_LENGTH_MAX-1) + "}"
        self.regExp = QRegExp(rx)
        v = RegExpValidator(editor, self.regExp)
        editor.setValidator(v)
        return editor


    def setEditorData(self, editor, index):
        editor.setAutoFillBackground(True)
        value = from_qvariant(index.model().data(index, Qt.DisplayRole), to_text_string)
        self.old_plabel = str(value)
        editor.setText(value)

#-------------------------------------------------------------------------------
# Line edit delegate for values in self.tableViewCoals
#-------------------------------------------------------------------------------

class ValueDelegate(QItemDelegate):
    def __init__(self, parent=None):
        super(ValueDelegate, self).__init__(parent)
        self.parent = parent


    def createEditor(self, parent, option, index):
        editor = QLineEdit(parent)
        if index.column() == 1 or index.column() == 2:
            v = DoubleValidator(editor, min=0.)
            v.setExclusiveMin(True)
        else:
            v = DoubleValidator(editor)
        editor.setValidator(v)
        return editor


    def setEditorData(self, editor, index):
        editor.setAutoFillBackground(True)
        value = from_qvariant(index.model().data(index, Qt.DisplayRole), to_text_string)
        editor.setText(value)


    def setModelData(self, editor, model, index):
        if not editor.isModified():
            return
        if editor.validator().state == QValidator.Acceptable:
            value = from_qvariant(editor.text(), float)
            for idx in self.parent.selectionModel().selectedIndexes():
                if idx.column() == index.column():
                    model.setData(idx, value, Qt.DisplayRole)

#-------------------------------------------------------------------------------
# StandarItemModel for Coals
#-------------------------------------------------------------------------------

class StandardItemModelCoals(QStandardItemModel):
    def __init__(self, case, model):
        """
        """
        QStandardItemModel.__init__(self)
        self.headers = [self.tr("Name"),
                        self.tr("Limit temperature\nof fouling (deg C)"),
                        self.tr("Ash critical\nviscosity (Pa.s)"),
                        self.tr("Coefficient 1"),
                        self.tr("Coefficient 2")]

        self.kwords = ["", "TPRENC", "VISREF", "ENC1" ,"ENC2"]

        self.setColumnCount(len(self.headers))
        self.dataCoals = []

        self.case = case
        self.model = model

        self.coalModel = CoalCombustionModel(self.case)
        CoalsNumber = self.coalModel.getCoalNumber()

        for icoal in range(CoalsNumber):
            line = []
            line.append(self.coalModel.getFuelLabel(icoal+1))
            line.append(self.model.getThresholdTemperatureOfFouling(icoal+1))
            line.append(self.model.getCriticalViscosityOfFouling(icoal+1))
            line.append(self.model.getCoef1OfFouling(icoal+1))
            line.append(self.model.getCoef2OfFouling(icoal+1))
            self.dataCoals.append(line)
            row = self.rowCount()
            self.setRowCount(row+1)


    def data(self, index, role):
        if not index.isValid():
            return None

        # ToolTips
        if role == Qt.ToolTipRole:
            return self.tr("Code_Saturne key word: " + self.kwords[index.column()])

        # Display
        if role == Qt.DisplayRole:
            return self.dataCoals[index.row()][index.column()]

        return None


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled
        if index.column() == 0:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable
        else:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable


    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return self.headers[section]
        return None


    def setData(self, index, value, role):
        row = index.row()
        col = index.column()

        val = from_qvariant(value, float)
        self.dataCoals[row][col] = val

        if col == 1:
            self.model.setThresholdTemperatureOfFouling(row+1, val)
        elif col == 2:
            self.model.setCriticalViscosityOfFouling(row+1, val)
        elif col == 3:
            self.model.setCoef1OfFouling(row+1, val)
        elif col == 4:
            self.model.setCoef2OfFouling(row+1, val)

        self.dataChanged.emit(index, index)
        return True


#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------


class LagrangianView(QWidget, Ui_LagrangianForm):
    """
    """
    def __init__(self, parent, case):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_LagrangianForm.__init__(self)
        self.setupUi(self)

        self.case = case
        self.case.undoStopGlobal()
        self.model = LagrangianModel(self.case)

        # Combo model
        self.modelIPHYLA = ComboModel(self.comboBoxIPHYLA,2,1)
        self.modelIPHYLA.addItem(self.tr("No model"), 'off')
        self.modelIPHYLA.addItem(self.tr("Heat transfer and evaporation"), 'thermal')
        if CoalCombustionModel(self.case).getCoalCombustionModel("only") != 'off':
            self.modelIPHYLA.addItem(self.tr("Pulverised coal model"), 'coal')

        # Connections
        self.checkBoxISUILA.clicked.connect(self.slotISUILA)
        self.checkBoxISTTIO.clicked.connect(self.slotISTTIO)
        self.checkBoxIDEPST.clicked.connect(self.slotIDEPST)
        self.comboBoxIPHYLA.activated[str].connect(self.slotIPHYLA)
        self.checkBoxITPVAR.clicked.connect(self.slotITPVAR)
        self.checkBoxIMPVAR.clicked.connect(self.slotIMPVAR)
        self.checkBoxIENCRA.clicked.connect(self.slotIENCRA)
        #
        self.lineEditNSTITS.textChanged[str].connect(self.slotNSTITS)
        self.checkBoxLTSDYN.clicked.connect(self.slotLTSDYN)
        self.checkBoxLTSMAS.clicked.connect(self.slotLTSMAS)
        self.checkBoxLTSTHE.clicked.connect(self.slotLTSTHE)
        self.toolButtonAdvanced.clicked.connect(self.slotAdvancedOptions)

        # Validators
        validatorNSTITS = IntValidator(self.lineEditNSTITS)

        self.lineEditNSTITS.setValidator(validatorNSTITS)

        # initialize Widgets
        model = self.model.getLagrangianModel()

        status = self.model.getRestart()
        if status == "on":
            self.checkBoxISUILA.setChecked(True)
        else:
            self.checkBoxISUILA.setChecked(False)

        status = self.model.getCarrierFlowStationary()
        if status == "on":
            self.checkBoxISTTIO.setChecked(True)
        else:
            self.checkBoxISTTIO.setChecked(False)

        status = self.model.getDepositionSubmodel()
        if status == "on":
            self.checkBoxIDEPST.setChecked(True)
        else:
            self.checkBoxIDEPST.setChecked(False)

        if ( model == "frozen" ):
            self.labelISTTIO.setDisabled(True)
            self.checkBoxISTTIO.setChecked(True)
            self.checkBoxISTTIO.setDisabled(True)

        self.groupBox2way.hide()

        self.labelISTTIO.setDisabled(False)
        self.checkBoxISTTIO.setDisabled(False)

        if model == "one_way":
            pass

        elif model == "two_way":
            self.groupBox2way.show()

            start_it = self.model.get2WayCouplingStartIteration()
            self.lineEditNSTITS.setText(str(start_it))

            status = self.model.get2WayCouplingDynamic()
            if status == "on":
                self.checkBoxLTSDYN.setChecked(True)
            else:
                self.checkBoxLTSDYN.setChecked(False)

            status = self.model.get2WayCouplingMass()
            if status == "on":
                self.checkBoxLTSMAS.setChecked(True)
            else:
                self.checkBoxLTSMAS.setChecked(False)

            status = self.model.get2WayCouplingTemperature()
            if status == "on":
                self.checkBoxLTSTHE.setChecked(True)
            else:
                self.checkBoxLTSTHE.setChecked(False)

        elif model == "frozen":
            self.labelISTTIO.setDisabled(True)
            self.checkBoxISTTIO.setDisabled(True)

        part_model = self.model.getParticlesModel()
        self.modelIPHYLA.setItem(str_model=part_model)
        self.slotIPHYLA(self.modelIPHYLA.dicoM2V[part_model])

        self.case.undoStartGlobal()


    @pyqtSlot()
    def slotISUILA(self):
        """
        Input ISUILA.
        """
        if self.checkBoxISUILA.isChecked():
            self.model.setRestart("on")
        else:
            self.model.setRestart("off")


    @pyqtSlot()
    def slotISTTIO(self):
        """
        Input ISTTIO.
        """
        if self.checkBoxISTTIO.isChecked():
            self.model.setCarrierFlowStationary("on")
        else:
            self.model.setCarrierFlowStationary("off")


    @pyqtSlot()
    def slotIDEPST(self):
        """
        Input IDEPST.
        """
        if self.checkBoxIDEPST.isChecked():
            self.model.setDepositionSubmodel("on")
        else:
            self.model.setDepositionSubmodel("off")


    @pyqtSlot(str)
    def slotIPHYLA(self, text):
        """
        Input IPHYLA.
        """
        value = self.modelIPHYLA.dicoV2M[str(text)]
        self.model.setParticlesModel(value)

        self.frameModel1.hide()
        self.frameModel2.hide()

        # No model
        if value == "off":
            pass

        # Equations on temperature, diameter and mass
        elif value == "thermal":

            self.frameModel1.show()

            status = self.model.getHeating()
            if status == "on":
                self.checkBoxITPVAR.setChecked(True)

            else:
                self.checkBoxITPVAR.setChecked(False)

            status = self.model.getEvaporation()
            if status == "on":
                self.checkBoxIMPVAR.setChecked(True)
            else:
                self.checkBoxIMPVAR.setChecked(False)

        # Pulverised coal model
        elif value == "coal":

            self.frameModel2.show()
            self.tableViewCoals.show()
            status = self.model.getCoalFouling()
            if status == "on":
                self.checkBoxIENCRA.setChecked(True)
            else:
                self.checkBoxIENCRA.setChecked(False)
            self.slotIENCRA()


    @pyqtSlot()
    def slotITPVAR(self):
        """
        Input ITPVAR.
        """
        if self.checkBoxITPVAR.isChecked():
            self.model.setHeating("on")
        else:
            self.model.setHeating("off")


    @pyqtSlot()
    def slotIMPVAR(self):
        """
        Input IMPVAR.
        """
        if self.checkBoxIMPVAR.isChecked():
            self.model.setEvaporation("on")
        else:
            self.model.setEvaporation("off")


    @pyqtSlot()
    def slotIENCRA(self):
        """
        Input IENCRA.
        """
        if self.checkBoxIENCRA.isChecked():
            self.model.setCoalFouling("on")

            self.modelCoals = StandardItemModelCoals(self.case, self.model)
            self.tableViewCoals.setModel(self.modelCoals)
            delegateValue = ValueDelegate(self.tableViewCoals)
            delegateValue2 = ValueDelegate(self.tableViewCoals)
            delegateLabel = LabelDelegate(self.tableViewCoals)
            self.tableViewCoals.setItemDelegateForColumn(0, delegateLabel)
            self.tableViewCoals.setItemDelegateForColumn(1, delegateValue)
            self.tableViewCoals.setItemDelegateForColumn(2, delegateValue)
            self.tableViewCoals.setItemDelegateForColumn(3, delegateValue2)
            self.tableViewCoals.setItemDelegateForColumn(4, delegateValue2)
            self.tableViewCoals.show()
            if QT_API == "PYQT4":
                self.tableViewCoals.verticalHeader().setResizeMode(QHeaderView.ResizeToContents)
                self.tableViewCoals.horizontalHeader().setResizeMode(QHeaderView.ResizeToContents)
                self.tableViewCoals.horizontalHeader().setResizeMode(0, QHeaderView.Stretch)
            elif QT_API == "PYQT5":
                self.tableViewCoals.verticalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
                self.tableViewCoals.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
                self.tableViewCoals.horizontalHeader().setSectionResizeMode(0, QHeaderView.Stretch)
        else:
            self.model.setCoalFouling("off")
            if hasattr(self, "modelCoals"):
                del self.modelCoals
            self.tableViewCoals.hide()


    @pyqtSlot(str)
    def slotNSTITS(self, text):
        """
        Input NSTITS.
        """
        if self.lineEditNSTITS.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, int)
            self.model.set2WayCouplingStartIteration(value)


    @pyqtSlot()
    def slotLTSDYN(self):
        """
        Input LTSDYN.
        """
        if self.checkBoxLTSDYN.isChecked():
            self.model.set2WayCouplingDynamic("on")
        else:
            self.model.set2WayCouplingDynamic("off")


    @pyqtSlot()
    def slotLTSMAS(self):
        """
        Input LTSMAS.
        """
        if self.checkBoxLTSMAS.isChecked():
            self.model.set2WayCouplingMass("on")
        else:
            self.model.set2WayCouplingMass("off")


    @pyqtSlot()
    def slotLTSTHE(self):
        """
        Input LTSTHE.
        """
        if self.checkBoxLTSTHE.isChecked():
            self.model.set2WayCouplingTemperature("on")
        else:
            self.model.set2WayCouplingTemperature("off")


    @pyqtSlot()
    def slotAdvancedOptions(self):
        """
        Ask one popup for advanced specifications
        """
        default = {}
        default['scheme_order']                        = self.model.getSchemeOrder()
        default['turbulent_dispertion']                = self.model.getTurbulentDispersion()
        default['fluid_particles_turbulent_diffusion'] = self.model.getTurbulentDiffusion()
        default['complete_model_iteration']            = self.model.getCompleteModelStartIteration()
        default['complete_model_direction']            = self.model.getCompleteModelDirection()

        dialog = LagrangianAdvancedOptionsDialogView(self, self.case, default)
        if dialog.exec_():
            result = dialog.get_result()
            self.model.setSchemeOrder(int(result['scheme_order']))
            self.model.setTurbulentDispersion(result['turbulent_dispertion'])
            self.model.setTurbulentDiffusion(result['fluid_particles_turbulent_diffusion'])
            self.model.setCompleteModelStartIteration(result['complete_model_iteration'])
            self.model.setCompleteModelDirection(int(result['complete_model_direction']))


    def tr(self, text):
        """
        Translation
        """
        return text

#-------------------------------------------------------------------------------
# Testing part
#-------------------------------------------------------------------------------

if __name__ == "__main__":
    pass

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
