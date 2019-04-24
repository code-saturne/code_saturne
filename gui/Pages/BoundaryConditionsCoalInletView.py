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
- BoundaryConditionsCoalInletView
- ValueDelegate
- StandardItemModelCoal
- StandardItemModelCoalMass
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
from code_saturne.Base.QtPage import DoubleValidator, ComboModel
from code_saturne.Base.QtPage import to_qvariant, from_qvariant, to_text_string

from code_saturne.Pages.BoundaryConditionsCoalInletForm import Ui_BoundaryConditionsCoalInletForm
import code_saturne.model.CoalCombustionModel as CoalCombustion

from code_saturne.model.LocalizationModel import LocalizationModel, Zone
from code_saturne.model.Boundary import Boundary
from code_saturne.Pages.QMegEditorView import QMegEditorView
from code_saturne.model.NotebookModel import NotebookModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("BoundaryConditionsCoalInletView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Line edit delegate with a Double validator (positive value)
#-------------------------------------------------------------------------------

class ValueDelegate(QItemDelegate):
    def __init__(self, parent=None):
        super(ValueDelegate, self).__init__(parent)
        self.parent = parent

    def createEditor(self, parent, option, index):
        editor = QLineEdit(parent)
        validator = DoubleValidator(editor, min=0.)
        editor.setValidator(validator)
        return editor

    def setEditorData(self, editor, index):
        editor.setAutoFillBackground(True)
        value = from_qvariant(index.model().data(index, Qt.DisplayRole), to_text_string)
        editor.setText(value)

    def setModelData(self, editor, model, index):
        if editor.validator().state == QValidator.Acceptable:
            value = from_qvariant(editor.text(), float)
            model.setData(index, to_qvariant(value), Qt.DisplayRole)

#-------------------------------------------------------------------------------
# StandarItemModel class to display Coals in a QTableView
#-------------------------------------------------------------------------------

class StandardItemModelCoal(QStandardItemModel):
    def __init__(self, case):
        QStandardItemModel.__init__(self)
        self.headers = [self.tr("Coal number"),
                        self.tr("Flow (kg/s)"),
                        self.tr("Temperature \n(K)")]
        self.setColumnCount(len(self.headers))
        self.dataCoal = []
        self.case = case


    def setBoundaryFromLabel(self, label):
        self.modelBoundary = Boundary('coal_inlet', label, self.case)


    def data(self, index, role):
        if not index.isValid():
            return to_qvariant()
        if role == Qt.DisplayRole:
            return to_qvariant(self.dataCoal[index.row()][index.column()])
        return to_qvariant()


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled
        elif index.column() in [1,2]:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable
        else:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable


    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return to_qvariant(self.headers[section])
        return to_qvariant()


    def setData(self, index, value, role):
        row = index.row()
        col = index.column()
        if not hasattr(self, "modelBoundary"):
            log.debug("ERROR in setData (StandardItemModelCoal) : no Boundary model defined")
            return
        v = from_qvariant(value, float)
        self.dataCoal[row][col] = v
        if col == 1:
            self.modelBoundary.setCoalFlow(v, row)
        elif col == 2:
            self.modelBoundary.setCoalTemperature(v, row)
        self.dataChanged.emit(index, index)
        return True


    def insertItem(self, nameCoal, valCoal, valCoalTemp):
        line = [nameCoal, valCoal, valCoalTemp]
        self.dataCoal.append(line)
        row = self.rowCount()
        self.setRowCount(row+1)


    def deleteAll(self):
        self.dataCoal = []
        self.setRowCount(0)

#-------------------------------------------------------------------------------
# StandarItemModel class to display Coal masses in a QTableView
#-------------------------------------------------------------------------------

class StandardItemModelCoalMass(QStandardItemModel):

    def __init__(self, case, coalNumber, coalClassesNumber):
        QStandardItemModel.__init__(self)
        self.case = case
        self.coalNumber = coalNumber
        self.coalClassesNumber = coalClassesNumber


    def setRatio(self, ratio):
        cols = len(ratio)
        if type(ratio[0]) == type([]):
            rows = max([len(c) for c in ratio])
        else:
            rows = 1
        self.setColumnCount(cols)
        self.setRowCount(rows)
        self.ratio = ratio


    def setBoundaryFromLabel(self, label):
        log.debug("setBoundaryFromLabel")
        self.modelBoundary = Boundary('coal_inlet', label, self.case)


    def data(self, index, role):
        if not index.isValid():
            return to_qvariant()
        if role == Qt.DisplayRole:
            classe = index.row()
            coal   = index.column()
            if classe < self.coalClassesNumber[coal]:
                try:
                    return to_qvariant(self.ratio[coal][classe])
                except:
                    log.debug("ERROR no data for self.ratio[%i][%i] "%(coal, classe))
        return to_qvariant()


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled
        elif index.row() >= self.coalClassesNumber[index.column()]:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable
        else:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable


    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return to_qvariant("Coal" + " " + str(section+1))
        if orientation == Qt.Vertical and role == Qt.DisplayRole:
            return to_qvariant("Class" + " " + str(section+1))
        return to_qvariant()


    def setData(self, index, value, role):
        if not hasattr(self, "modelBoundary"):
            log.debug("ERROR in setData (StandardItemModelCoalMass): no Boundary model defined")
            return
        classe = index.row()
        coal   = index.column()
        if not value:
            return False
        v = from_qvariant(value, float)
        self.ratio[coal][classe] = v
        log.debug("setData v = %f "%v)

        lst = self.modelBoundary.getCoalRatios(coal)
        lastValue = 0
        for iclasse in range(0, self.coalClassesNumber[coal]-1):
            lastValue += self.ratio[coal][iclasse]

        if lastValue < 100.+ 1e-6 :
            lst[classe] = self.ratio[coal][classe]
            lastValue = 100 - lastValue
            self.ratio[coal][self.coalClassesNumber[coal]-1] = lastValue
            lst[self.coalClassesNumber[coal]-1] = lastValue
            self.modelBoundary.setCoalRatios(coal, lst)
        else :
            self.ratio[coal][classe] = lst[classe]

        self.dataChanged.emit(index, index)
        return True


    def deleteAll(self):
        self.ratio = []
        self.setRowCount(0)

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class BoundaryConditionsCoalInletView(QWidget, Ui_BoundaryConditionsCoalInletForm):
    def __init__(self, parent):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_BoundaryConditionsCoalInletForm.__init__(self)
        self.setupUi(self)


    def setup(self, case):
        """
        Setup the widget
        """
        self.case = case
        self.__boundary = None

        self.case.undoStopGlobal()
        self.notebook = NotebookModel(self.case)

        # Connections
        self.comboBoxTypeInlet.activated[str].connect(self.__slotInletType)
        self.comboBoxVelocity.activated[str].connect(self.__slotChoiceVelocity)
        self.lineEditVelocity.textChanged[str].connect(self.__slotVelocityValue)
        self.lineEditTemperature.textChanged[str].connect(self.__slotTemperature)
        self.spinBoxOxydantNumber.valueChanged[int].connect(self.__slotOxydantNumber)

        self.comboBoxDirection.activated[str].connect(self.__slotChoiceDirection)
        self.lineEditDirectionX.textChanged[str].connect(self.__slotDirX)
        self.lineEditDirectionY.textChanged[str].connect(self.__slotDirY)
        self.lineEditDirectionZ.textChanged[str].connect(self.__slotDirZ)

        # Combo models
        self.modelTypeInlet = ComboModel(self.comboBoxTypeInlet, 2, 1)
        self.modelTypeInlet.addItem(self.tr("only oxydant"), 'oxydantFlow')
        self.modelTypeInlet.addItem(self.tr("oxydant and coal"), 'coalFlow')

        self.modelVelocity = ComboModel(self.comboBoxVelocity, 4, 1)
        self.modelVelocity.addItem(self.tr("norm"), 'norm')
        self.modelVelocity.addItem(self.tr("mass flow rate"), 'flow1')
        self.modelVelocity.addItem(self.tr("norm (user law)"), 'norm_formula')
        self.modelVelocity.addItem(self.tr("mass flow rate (user law)"), 'flow1_formula')

        self.modelDirection = ComboModel(self.comboBoxDirection, 3, 1)
        self.modelDirection.addItem(self.tr("normal to the inlet"), 'normal')
        self.modelDirection.addItem(self.tr("specified coordinates"), 'coordinates')
        self.modelDirection.addItem(self.tr("user profile"), 'formula')

        # Validators
        validatorVelocity = DoubleValidator(self.lineEditVelocity)
        validatorX = DoubleValidator(self.lineEditDirectionX)
        validatorY = DoubleValidator(self.lineEditDirectionY)
        validatorZ = DoubleValidator(self.lineEditDirectionZ)
        validatorTemp = DoubleValidator(self.lineEditTemperature, min=0.)

        # Apply validators
        self.lineEditVelocity.setValidator(validatorVelocity)
        self.lineEditDirectionX.setValidator(validatorX)
        self.lineEditDirectionY.setValidator(validatorY)
        self.lineEditDirectionZ.setValidator(validatorZ)
        self.lineEditTemperature.setValidator(validatorTemp)

        self.pushButtonVelocityFormula.clicked.connect(self.__slotVelocityFormula)
        self.pushButtonDirectionFormula.clicked.connect(self.__slotDirectionFormula)

        # Useful information about coals, classes, and ratios

        mdl =  CoalCombustion.CoalCombustionModel(self.case)
        if mdl.getCoalCombustionModel() != "off":
            self.__coalNumber = mdl.getCoalNumber()
            self.__coalClassesNumber = []
            for coal in range(0, self.__coalNumber):
                self.__coalClassesNumber.append(mdl.getClassNumber(str(coal+1)))
            self.__maxOxydantNumber = mdl.getOxidantNumber()
        else:
            self.__coalNumber = 0
            self.__coalClassesNumber = [0]
            self.__maxOxydantNumber = 1

        self.__ratio = self.__coalNumber*[0]
        for i in range(0, self.__coalNumber):
            self.__ratio[i] = self.__coalClassesNumber[i]*[0]

        # Coal table

        self.__modelCoal = StandardItemModelCoal(self.case)
        self.tableViewCoal.setModel(self.__modelCoal)
        delegateValue = ValueDelegate(self.tableViewCoal)
        self.tableViewCoal.setItemDelegateForColumn(1, delegateValue)
        self.tableViewCoal.setItemDelegateForColumn(2, delegateValue)

        # Coal mass ratio table

        self.__modelCoalMass = StandardItemModelCoalMass(self.case,
                                                         self.__coalNumber,
                                                         self.__coalClassesNumber)
        self.tableViewCoalMass.setModel(self.__modelCoalMass)

        delegateValueMass = ValueDelegate(self.tableViewCoalMass)
        for c in range(self.__modelCoalMass.columnCount()):
            self.tableViewCoalMass.setItemDelegateForColumn(c, delegateValueMass)

        self.case.undoStartGlobal()


    def showWidget(self, b):
        """
        Show the widget
        """
        label = b.getLabel()
        self.__boundary = Boundary('coal_inlet', label, self.case)

        # Initialize velocity
        choice = self.__boundary.getVelocityChoice()
        self.modelVelocity.setItem(str_model=choice)
        self.__updateLabel()

        if choice[-7:] == "formula":
            self.pushButtonVelocityFormula.setEnabled(True)
            self.lineEditVelocity.setEnabled(False)
        else:
            self.pushButtonVelocityFormula.setEnabled(False)
            self.lineEditVelocity.setEnabled(True)
            v = self.__boundary.getVelocity()
            self.lineEditVelocity.setText(str(v))

        # Initialize oxydant and temperature
        self.spinBoxOxydantNumber.setMaximum(self.__maxOxydantNumber)
        o = self.__boundary.getOxydantNumber()
        self.spinBoxOxydantNumber.setValue(o)
        t = self.__boundary.getOxydantTemperature()
        self.lineEditTemperature.setText(str(t))

        # Initialize direction
        choice = self.__boundary.getDirectionChoice()
        self.modelDirection.setItem(str_model=choice)
        text = self.modelDirection.dicoM2V[choice]
        if choice == "formula":
            self.pushButtonDirectionFormula.setEnabled(True)
            self.frameDirectionCoordinates.hide()
        elif choice == "coordinates":
            self.pushButtonDirectionFormula.setEnabled(False)
            self.frameDirectionCoordinates.show()
            v = self.__boundary.getDirection('direction_x')
            self.lineEditDirectionX.setText(str(v))
            v = self.__boundary.getDirection('direction_y')
            self.lineEditDirectionY.setText(str(v))
            v = self.__boundary.getDirection('direction_z')
            self.lineEditDirectionZ.setText(str(v))
        elif choice == "normal":
            self.pushButtonDirectionFormula.setEnabled(False)
            self.frameDirectionCoordinates.hide()

        log.debug("showWidget:inlet type: %s " % self.__boundary.getInletType())
        if self.__boundary.getInletType() == "coalFlow":
            self.modelTypeInlet.setItem(str_model="coalFlow")
            self.groupBoxCoal.show()
            self.groupBoxCoalMass.show()
            self.__updateTables()
            self.__boundary.setInletType("coalFlow")
        else:
            self.__boundary.setInletType("oxydantFlow")
            self.modelTypeInlet.setItem(str_model="oxydantFlow")
            self.groupBoxCoal.hide()
            self.groupBoxCoalMass.hide()

        self.show()


    def hideWidget(self):
        """
        Hide all
        """
        self.hide()


    def __updateTables(self):
        """
        Insert rows in the two QTableView.
        """
        # clean the QTableView
        self.__modelCoal.deleteAll()
        self.__modelCoalMass.deleteAll()

        label = self.__boundary.getLabel()
        self.__modelCoalMass.setBoundaryFromLabel(label)
        self.__modelCoal.setBoundaryFromLabel(label)

        # fill the flow and temperature of the coal
        for coal in range(0, self.__coalNumber):
            self.__modelCoal.insertItem(self.tr("Coal ") + " " + str(coal+1),
                                        self.__boundary.getCoalFlow(coal),
                                        self.__boundary.getCoalTemperature(coal))

        # fill the ratio of mass for each class for each coal
        for coal in range(0, self.__coalNumber) :
            lastValue = 0.
            for coalClass in range(0, self.__coalClassesNumber[coal]-1):
                lst = self.__boundary.getCoalRatios(coal)
                lastValue += lst[coalClass]
                self.__ratio[coal][coalClass] = lst[coalClass]

            # last class is computed in order to assure that sum is egal to 100%
            coalClass = self.__coalClassesNumber[coal]-1
            lastValue = 100 - lastValue
            self.__ratio[coal][coalClass] = lastValue

        self.__modelCoalMass.setRatio(self.__ratio)


    @pyqtSlot(str)
    def __slotChoiceVelocity(self, text):
        """
        Private slot.

        Input the velocity boundary type choice (norm, ).

        @type text: C{QString}
        @param text: velocity boundary type choice.
        """
        c = self.modelVelocity.dicoV2M[str(text)]
        log.debug("slotChoiceVelocity: %s " % c)
        self.__boundary.setVelocityChoice(c)

        if c[-7:] == "formula":
            self.pushButtonVelocityFormula.setEnabled(True)
            exp = self.__boundary.getVelocity()
            if exp:
                self.pushButtonVelocityFormula.setStyleSheet("background-color: green")
                self.pushButtonVelocityFormula.setToolTip(exp)
            else:
                self.pushButtonVelocityFormula.setStyleSheet("background-color: red")
            self.lineEditVelocity.setEnabled(False)
            self.lineEditVelocity.setText("")
        else:
            self.pushButtonVelocityFormula.setEnabled(False)
            self.pushButtonVelocityFormula.setStyleSheet("background-color: None")
            self.lineEditVelocity.setEnabled(True)
            v = self.__boundary.getVelocity()
            self.lineEditVelocity.setText(str(v))

        self.__updateLabel()


    def __updateLabel(self):
        """
        Update the unit for the velocity specification.
        """
        c = self.__boundary.getVelocityChoice()
        if c in ('norm', 'norm_formula'):
            self.labelUnitVelocity.setText(str('m/s'))
        elif c in ('flow1', 'flow1_formula'):
            self.labelUnitVelocity.setText(str('kg/s'))
        elif c in ('flow2', 'flow2_formula'):
            self.labelUnitVelocity.setText(str('m<sup>3</sup>/s'))


    @pyqtSlot(str)
    def __slotVelocityValue(self, text):
        """
        Private slot.

        New value associated to the velocity boundary type.

        @type text: C{QString}
        @param text: value
        """
        if self.sender().validator().state == QValidator.Acceptable:
            v = from_qvariant(text, float)
            self.__boundary.setVelocity(v)


    @pyqtSlot()
    def __slotVelocityFormula(self):
        """
        """
        exp = self.__boundary.getVelocity()
        c = self.__boundary.getVelocityChoice()
        if c == 'norm_formula':
            exa = "u_norm = 1.0;"
            req = [('u_norm', 'Norm of the velocity')]
        elif c == 'flow1_formula':
            exa = "q_m = 1.0;"
            req = [('q_m', 'mass flow rate')]
        elif c == 'flow2_formula':
            exa = "q_v = 1.0;"
            req = [('q_v', 'volumic flow rate')]

        sym = [('x', "X face's gravity center"),
               ('y', "Y face's gravity center"),
               ('z', "Z face's gravity center"),
               ('dt', 'time step'),
               ('t', 'current time'),
               ('iter', 'number of iteration')]

        for (nme, val) in self.notebook.getNotebookList():
            sym.append((nme, 'value (notebook) = ' + str(val)))

        dialog = QMegEditorView(parent = self,
                                function_type = 'bnd',
                                zone_name     = self.__boundary._label,
                                variable_name = 'velocity',
                                expression    = exp,
                                required      = req,
                                symbols       = sym,
                                condition     = c,
                                examples      = exa)

        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormulaVelocity -> %s" % str(result))
            self.__boundary.setVelocity(str(result))
            self.pushButtonVelocityFormula.setToolTip(result)
            self.pushButtonVelocityFormula.setStyleSheet("background-color: green")


    @pyqtSlot(str)
    def __slotChoiceDirection(self, text):
        """
        Input the direction type choice.
        """
        c = self.modelDirection.dicoV2M[str(text)]
        log.debug("slotChoiceVelocity: %s " % c)
        self.__boundary.setDirectionChoice(c)

        if c == "formula":
            self.pushButtonDirectionFormula.setEnabled(True)
            exp = self.__boundary.getDirection('direction_formula')
            if exp:
                self.pushButtonDirectionFormula.setStyleSheet("background-color: green")
                self.pushButtonDirectionFormula.setToolTip(exp)
            else:
                self.pushButtonDirectionFormula.setStyleSheet("background-color: red")
            self.frameDirectionCoordinates.hide()
        elif c == "coordinates":
            self.pushButtonDirectionFormula.setEnabled(False)
            self.pushButtonDirectionFormula.setStyleSheet("background-color: None")
            self.frameDirectionCoordinates.show()
            v = self.__boundary.getDirection('direction_x')
            self.lineEditDirectionX.setText(str(v))
            v = self.__boundary.getDirection('direction_y')
            self.lineEditDirectionY.setText(str(v))
            v = self.__boundary.getDirection('direction_z')
            self.lineEditDirectionZ.setText(str(v))
        elif c == "normal":
            self.pushButtonDirectionFormula.setEnabled(False)
            self.pushButtonDirectionFormula.setStyleSheet("background-color: None")
            self.frameDirectionCoordinates.hide()


    @pyqtSlot(str)
    def __slotDirX(self, text):
        """
        INPUT value into direction of inlet flow
        """
        if self.sender().validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.__boundary.setDirection('direction_x', value)


    @pyqtSlot(str)
    def __slotDirY(self, text):
        """
        INPUT value into direction of inlet flow
        """
        if self.sender().validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.__boundary.setDirection('direction_y', value)


    @pyqtSlot(str)
    def __slotDirZ(self, text):
        """
        INPUT value into direction of inlet flow
        """
        if self.sender().validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.__boundary.setDirection('direction_z', value)


    @pyqtSlot()
    def __slotDirectionFormula(self):
        """
        """
        exp = self.__boundary.getDirection('direction_formula')

        req = [('dir_x', 'Direction of the flow along X'),
               ('dir_y', 'Direction of the flow along Y'),
               ('dir_z', 'Direction of the flow along Z')]

        exa = "dir_x = 3.0;\ndir_y = 1.0;\ndir_z = 0.0;\n"

        sym = [('x', "X face's gravity center"),
               ('y', "Y face's gravity center"),
               ('z', "Z face's gravity center"),
               ('dt', 'time step'),
               ('t', 'current time'),
               ('iter', 'number of iteration')]

        for (nme, val) in self.notebook.getNotebookList():
            sym.append((nme, 'value (notebook) = ' + str(val)))

        dialog = QMegEditorView(parent = self,
                                function_type = 'bnd',
                                zone_name     = self.__boundary._label,
                                variable_name = 'direction',
                                expression    = exp,
                                required      = req,
                                symbols       = sym,
                                condition     = 'formula',
                                examples      = exa)

        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormulaDirection -> %s" % str(result))
            self.__boundary.setDirection('direction_formula', str(result))
            self.pushButtonDirectionFormula.setToolTip(result)
            self.pushButtonDirectionFormula.setStyleSheet("background-color: green")



    @pyqtSlot(str)
    def __slotInletType(self, text):
        """
        INPUT inlet type : 'oxydant' or 'oxydant + coal'
        """
        value = self.modelTypeInlet.dicoV2M[str(text)]
        log.debug("__slotInletType value = %s " % value)

        self.__boundary.setInletType(value)

        if value == 'oxydantFlow':
            self.groupBoxCoal.hide()
            self.groupBoxCoalMass.hide()
        else:
            self.groupBoxCoal.show()
            self.groupBoxCoalMass.show()
            self.__updateTables()


    @pyqtSlot(str)
    def __slotTemperature(self, text):
        if self.sender().validator().state == QValidator.Acceptable:
            t = from_qvariant(text, float)
            self.__boundary.setOxydantTemperature(t)


    @pyqtSlot(int)
    def __slotOxydantNumber(self, i):
        self.__boundary.setOxydantNumber(i)


    def getCoalNumber(self):
        """
        Return the coal number
        """
        return self.__coalNumber


    def tr(self, text):
        """
        Translation
        """
        return text

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
