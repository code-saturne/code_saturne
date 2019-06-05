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
This module defines the Immersed boundaries view data management.

This module contains the following classes and function:
- SyrthesVerbosityDelegate
- ProjectionAxisDelegate
- SelectionCriteriaDelegate
- StandardItemModelSyrthes
- ConjugateHeatTransferView
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

from code_saturne.model.Common import LABEL_LENGTH_MAX, GuiParam
from code_saturne.Base.QtPage import IntValidator, DoubleValidator, RegExpValidator, ComboModel
from code_saturne.Base.QtPage import from_qvariant, to_text_string
from code_saturne.Pages.ImmersedBoundariesNeptune import Ui_ImmersedBoundariesNeptune
from code_saturne.model.ImmersedBoundariesModel import ImmersedBoundariesModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("ImmersedBoundariesViewNeptune")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# QLineEdit delegate to attach a label to the FSI object
#-------------------------------------------------------------------------------

class FSIObjectNameDelegate(QItemDelegate):

    def __init__(self, parent = None):
        super(FSIObjectNameDelegate, self).__init__(parent)
        self.parent = parent


    def createEditor(self, parent, option, index):
        editor = QLineEdit(parent)
        return editor


    def setEditorData(self, editor, index):
        self.value = from_qvariant(index.model().data(index, Qt.DisplayRole),
                                   to_text_string)
        editor.setText(self.value)


    def setModelData(self, editor, model, index):
        value = editor.text()

        if str(value) != "":
            model.setData(index, value, Qt.DisplayRole)


#-------------------------------------------------------------------------------
# QComboBox delegate for the FSI type : Set motion or computed from fluid forces
#-------------------------------------------------------------------------------

class FSITypeDelegate(QItemDelegate):
    """
    Use of a combobox to set the fsi interaction type
    """

    def __init__(self, parent, mdl):
        super(FSITypeDelegate, self).__init__(parent)
        self.parent  = parent
        self.mdl     = mdl


    def createEditor(self, parent, option, index):
        editor = QComboBox(parent)

        for itm in ["imposed", "computed"]:
            editor.addItem(itm)

        editor.installEventFilter(self)
        return editor


    def setEditorData(self, comboBox, index):
        row = index.row()
        col = index.column()
        string = index.model().dataFSI[row][col]
        comboBox.setEditText(string)


    def setModelData(self, comboBox, model, index):
        value = comboBox.currentText()
        model.setData(index, value, Qt.DisplayRole)

#-------------------------------------------------------------------------------
# QLineEdit delegate for validation of Syrthes verbosity or visualization
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# QComboBox delegate for Axis Projection in Conjugate Heat Transfer table
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# StandarItemModel class
#-------------------------------------------------------------------------------

class StandardItemModelFSI(QStandardItemModel):

    def __init__(self, model):
        """
        """
        QStandardItemModel.__init__(self)

        self.headers = [self.tr("Object name"),
                        self.tr("Interaction type")]
        self.tooltip = [self.tr("Name of solid object"),
                        self.tr("Type of motion interaction with the flow")]

        self.setColumnCount(len(self.headers))
        self.dataFSI = []
        self.__model = model


    def data(self, index, role):
        if not index.isValid():
            return None

        # Tooltips
        if role == Qt.ToolTipRole:
            return self.tooltip[index.column()]

        # Display
        if role == Qt.DisplayRole:
            return self.dataFSI[index.row()][index.column()]
        elif role == Qt.TextAlignmentRole:
            return Qt.AlignCenter

        return None


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled

        return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable


    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return self.headers[section]
        return None


    def setData(self, index, value, role):
        if not index.isValid():
            return

        row = index.row()
        col = index.column()

        self.dataFSI[row][col] = str(from_qvariant(value, to_text_string))

        num = row + 1
        self.__model.setObjectName(num, self.dataFSI[row][0])
        self.__model.setObjectInteraction(num, self.dataFSI[row][1])

#        self.dataChanged.emit(index, index)

        id1 = self.index(0, 0)
        id2 = self.index(self.rowCount(), 0)
        self.dataChanged.emit(id1, id2)
        return True


    def getData(self, index):
        row = index.row()
        return self.dataFSI[row]

    def addItem(self, object_name, interaction_type):
        """
        Add a row in the table.
        """
        self.dataFSI.append([object_name, interaction_type])
        row = self.rowCount()
        self.setRowCount(row+1)


    def deleteRow(self, row):
        """
        Delete the row in the model
        """
        del self.dataFSI[row]
        row = self.rowCount()
        self.setRowCount(row-1)

    def getItem(self, row):
        return self.dataFSI[row]

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class ImmersedBoundariesViewNeptune(QWidget, Ui_ImmersedBoundariesNeptune):
    """
    """
    def __init__(self, parent, case):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_ImmersedBoundariesNeptune.__init__(self)
        self.setupUi(self)

        self.case = case
        self.case.undoStopGlobal()

        self.ibm = ImmersedBoundariesModel(self.case)

        self.current_obj = None

        # Models
        self.modelFSI = StandardItemModelFSI(self.ibm)
        self.tableViewFSI.setModel(self.modelFSI)

        for obj in range(1,self.ibm.getNumberOfFSIObjects()+1):
            self.modelFSI.addItem(self.ibm.getObjectName(obj),
                                  self.ibm.getObjectInteraction(obj))

        if QT_API == "PYQT4":
            self.tableViewFSI.verticalHeader().setResizeMode(QHeaderView.ResizeToContents)
            self.tableViewFSI.horizontalHeader().setResizeMode(QHeaderView.ResizeToContents)
            self.tableViewFSI.horizontalHeader().setResizeMode(2, QHeaderView.Stretch)
        elif QT_API == "PYQT5":
            self.tableViewFSI.verticalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
            self.tableViewFSI.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
            self.tableViewFSI.horizontalHeader().setSectionResizeMode(2, QHeaderView.Stretch)

        self.modelFSI.dataChanged.connect(self.dataChanged)

        delegateObjectLabel  = FSIObjectNameDelegate(self.tableViewFSI)
        self.tableViewFSI.setItemDelegateForColumn(0, delegateObjectLabel)

        delegateObjectType   = FSITypeDelegate(self.tableViewFSI, self.ibm)
        self.tableViewFSI.setItemDelegateForColumn(1, delegateObjectType)

        self.checkBoxActivate.stateChanged.connect(self.slotCheckActivate)

        self.tableViewFSI.clicked[QModelIndex].connect(self.slotChangedSelection)


        for ind in ['Explicit', 'MEDCoupling']:
            eval('self.radioButton'+ind+'.toggled.connect(self.slotRadioButton)')

        # Connections
        self.pushButtonAddFSI.clicked.connect(self.slotAddFSI)
        self.pushButtonDeleteFSI.clicked.connect(self.slotDeleteFSI)

        self.pushButtonExplicit.clicked.connect(self.slotExplicitFormula)

        validatorDensity = DoubleValidator(self.lineEditObjDensity, min = 0.0)
        self.lineEditObjDensity.setValidator(validatorDensity)
        self.lineEditObjDensity.textChanged[str].connect(self.slotObjDensity)

        validatorStiffness = DoubleValidator(self.lineEditObjStiffness, min = 0.0)
        self.lineEditObjStiffness.setValidator(validatorStiffness)
        self.lineEditObjStiffness.textChanged[str].connect(self.slotObjStiffness)

        validatorDamping = DoubleValidator(self.lineEditObjDamping, min = 0.0)
        self.lineEditObjDamping.setValidator(validatorDamping)
        self.lineEditObjDamping.textChanged[str].connect(self.slotObjDamping)

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

        # Check for MEDCoupling presence
        import cs_config
        cfg = cs_config.config()
        self.has_medcoupling = cfg.libs['medcoupling'].have == 'yes'
        # deactivated for the moment
        self.has_medcoupling = False

        self.radioButtonMEDCoupling.setEnabled(self.has_medcoupling)
        if self.ibm.getMethod() == 'medcoupling' and self.has_medcoupling == False:
            self.setMethod('explicit')

        # Show/hide widgets on start
        if self.ibm.getOnOff() == 'off':
            self.groupBoxMethod.hide()
            self.groupBoxObjects.hide()
            self.groupBoxObjProperties.hide()
            self.groupBoxExplicit.hide()
            self.groupBoxMEDCoupling.hide()
        else:
            self.groupBoxMethod.show()
            self.groupBoxObjects.show()
            if self.ibm.getMethod() == 'explicit':
                self.groupBoxExplicit.show()
            else:
                self.groupBoxMEDCoupling.show()

        self.updatePageView()

        self.case.undoStartGlobal()


    @pyqtSlot("QModelIndex")
    def slotChangedSelection(self, index):
        """
        detect change in selection and update view
        """
        row = self.tableViewFSI.currentIndex().row()
        self.current_obj = row + 1
        self.updatePageView()


    @pyqtSlot(int)
    def slotCheckActivate(self, val):

        # Set the method state
        if val == 0:
            self.ibm.setOnOff('off')
        else:
            self.ibm.setOnOff('on')

        # Update the view if needed
        self.updatePageView()

    def dataChanged(self, topLeft, bottomRight):
        self.updatePageView()

    def updatePageView(self):

        if self.ibm.getOnOff() == 'off':
            self.checkBoxActivate.setChecked(False)
            self.groupBoxMethod.hide()
            self.groupBoxObjects.hide()
            self.groupBoxExplicit.hide()
            self.groupBoxMEDCoupling.hide()
            self.radioButtonExplicit.setChecked(False)
            self.radioButtonMEDCoupling.setChecked(False)

        else:
            self.checkBoxActivate.setChecked(True)
            self.groupBoxMethod.show()
            self.groupBoxObjects.show()

            # Which button to show for the solid definition
            if self.current_obj:
                if self.ibm.getMethod() == 'explicit':
                    self.groupBoxExplicit.show()
                    self.groupBoxMEDCoupling.hide()
                    self.radioButtonExplicit.setChecked(True)
                    self.radioButtonMEDCoupling.setChecked(False)
                elif self.ibm.getMethod() == 'medcoupling':
                    self.groupBoxExplicit.hide()
                    self.groupBoxMEDCoupling.show()
                    self.radioButtonExplicit.setChecked(False)
                    self.radioButtonMEDCoupling.setChecked(True)
                else:
                    self.radioButtonExplicit.setChecked(False)
                    self.radioButtonMEDCoupling.setChecked(False)
                    self.groupBoxExplicit.hide()
                    self.groupBoxMEDCoupling.hide()

                if self.ibm.getObjectInteraction(self.current_obj) == 'computed':
                    self.groupBoxObjProperties.show()
                    # Set correct values for each slot
                    self.lineEditObjDensity.setText(str(
                            self.ibm.getObjectDensity(self.current_obj)))

                    self.lineEditObjStiffness.setText(str(
                            self.ibm.getObjectStiffness(self.current_obj)))

                    self.lineEditObjDamping.setText(str(
                            self.ibm.getObjectDamping(self.current_obj)))

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
                else:
                    self.groupBoxObjProperties.hide()

            else:
                self.groupBoxExplicit.hide()
                self.groupBoxMEDCoupling.hide()
                self.groupBoxObjProperties.hide()


    @pyqtSlot()
    def slotRadioButton(self):

        for ind in ['Explicit', 'MEDCoupling']:

            radioButton = eval('self.radioButton'+ind)
            if radioButton.isChecked():
                self.ibm.setMethod(ind.lower())

        self.updatePageView()


    @pyqtSlot()
    def slotAddFSI(self):

        name        = '_'.join([self.ibm.defaultValues()['fsi_object_name'],
                                str(self.ibm.getNumberOfFSIObjects()+1)])
        interaction = self.ibm.defaultValues()['fsi_interaction']

        num = self.ibm.addFSIObject(name, interaction)
        self.modelFSI.addItem(name, interaction)


    @pyqtSlot()
    def slotDeleteFSI(self):
        row = self.tableViewFSI.currentIndex().row()
        log.debug("slotDeleteFSI -> %s" % (row,))
        if row == -1:
            title = self.tr("Warning")
            msg   = self.tr("You must select an existing object")
            QMessageBox.information(self, title, msg)
        else:
            self.modelFSI.deleteRow(row)
            self.ibm.deleteFSIObject(row+1)

    def tr(self, text):
        """
        Translation
        """
        return text

    @pyqtSlot()
    def slotExplicitFormula(self):
        """
        Explicit formula for variable porosity
        """

        objId = self.current_obj

        exp, req, sym = self.ibm.getFormulaPorosityComponents(objId)
        exa = ""

        name = self.ibm.getObjectName(objId)

        dialog = QMegEditorView(parent        = self,
                                function_type = 'var_poro',
                                zone_name     = name,
                                variable_name = 'porosity',
                                expression    = exp,
                                required      = req,
                                symbols       = sym,
                                known_fields  = [],
                                examples      = exa)

        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotExplicitFormula -> %s" % str(result))
            self.ibm.setExplicitFormula(str(objId), 'porosity', result)
            self.pushButtonExplicit.setStyleSheet("background-color: green")
            self.pushButtonExplicit.setToolTip(exp)


    @pyqtSlot(str)
    def slotObjDensity(self, text):
        num = self.tableViewFSI.currentIndex().row() + 1
        val = float(text)
        self.ibm.setObjectDensity(num, val)


    @pyqtSlot(str)
    def slotObjStiffness(self, text):
        num = self.tableViewFSI.currentIndex().row() + 1
        val = float(text)
        self.ibm.setObjectStiffness(num, val)


    @pyqtSlot(str)
    def slotObjDamping(self, text):
        num = self.tableViewFSI.currentIndex().row() + 1
        val = float(text)
        self.ibm.setObjectDamping(num, val)


    @pyqtSlot(str)
    def slotObjXinit(self, text):
        num = self.tableViewFSI.currentIndex().row() + 1
        val = float(text)
        self.ibm.setObjectInitPosition(num, xini=val)


    @pyqtSlot(str)
    def slotObjYinit(self, text):
        num = self.tableViewFSI.currentIndex().row() + 1
        val = float(text)
        self.ibm.setObjectInitPosition(num, yini=val)


    @pyqtSlot(str)
    def slotObjZinit(self, text):
        num = self.tableViewFSI.currentIndex().row() + 1
        val = float(text)
        self.ibm.setObjectInitPosition(num, zini=val)


    @pyqtSlot(str)
    def slotObjXeq(self, text):
        num = self.tableViewFSI.currentIndex().row() + 1
        val = float(text)
        self.ibm.setObjectEqPosition(num, xeq=val)


    @pyqtSlot(str)
    def slotObjYeq(self, text):
        num = self.tableViewFSI.currentIndex().row() + 1
        val = float(text)
        self.ibm.setObjectEqPosition(num, yeq=val)


    @pyqtSlot(str)
    def slotObjZeq(self, text):
        num = self.tableViewFSI.currentIndex().row() + 1
        val = float(text)
        self.ibm.setObjectEqPosition(num, zeq=val)


    @pyqtSlot(str)
    def slotObjVelXinit(self, text):
        num = self.tableViewFSI.currentIndex().row() + 1
        val = float(text)
        self.ibm.setObjectInitVel(num, vx=val)


    @pyqtSlot(str)
    def slotObjVelYinit(self, text):
        num = self.tableViewFSI.currentIndex().row() + 1
        val = float(text)
        self.ibm.setObjectInitVel(num, vy=val)


    @pyqtSlot(str)
    def slotObjVelZinit(self, text):
        num = self.tableViewFSI.currentIndex().row() + 1
        val = float(text)
        self.ibm.setObjectInitVel(num, vz=val)


    @pyqtSlot(str)
    def slotObjAccXinit(self, text):

        num = self.tableViewFSI.currentIndex().row() + 1
        val = float(text)
        self.ibm.setObjectInitAcc(num, ax=val)

    @pyqtSlot(str)
    def slotObjAccYinit(self, text):

        num = self.tableViewFSI.currentIndex().row() + 1
        val = float(text)
        self.ibm.setObjectInitAcc(num, ay=val)

    @pyqtSlot(str)
    def slotObjAccZinit(self, text):

        num = self.tableViewFSI.currentIndex().row() + 1
        val = float(text)
        self.ibm.setObjectInitAcc(num, az=val)


#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
