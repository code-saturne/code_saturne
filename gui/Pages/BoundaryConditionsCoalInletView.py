# -*- coding: iso-8859-1 -*-
#
#-------------------------------------------------------------------------------
#
#     This file is part of the Code_Saturne User Interface, element of the
#     Code_Saturne CFD tool.
#
#     Copyright (C) 1998-2009 EDF S.A., France
#
#     contact: saturne-support@edf.fr
#
#     The Code_Saturne User Interface is free software; you can redistribute it
#     and/or modify it under the terms of the GNU General Public License
#     as published by the Free Software Foundation; either version 2 of
#     the License, or (at your option) any later version.
#
#     The Code_Saturne User Interface is distributed in the hope that it will be
#     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
#     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with the Code_Saturne Kernel; if not, write to the
#     Free Software Foundation, Inc.,
#     51 Franklin St, Fifth Floor,
#     Boston, MA  02110-1301  USA
#
#-------------------------------------------------------------------------------

"""
This module contains the following classes:
-  BoundaryConditionsCoalInletView
- ValueDelegate
- StandardItemModelCoal
- StandardItemModelCoalMass
"""

#-------------------------------------------------------------------------------
# Standard modules
#-------------------------------------------------------------------------------

import string, logging

#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------

from PyQt4.QtCore import *
from PyQt4.QtGui  import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from Pages.BoundaryConditionsCoalInletForm import Ui_BoundaryConditionsCoalInletForm

from Base.Toolbox import GuiParam
from Base.QtPage import DoubleValidator, ComboModel
from Pages.LocalizationModel import LocalizationModel, Zone
from Pages.Boundary import Boundary

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
        #editor.installEventFilter(self)
        return editor

    def setEditorData(self, editor, index):
        value = index.model().data(index, Qt.DisplayRole).toString()
        editor.setText(value)

    def setModelData(self, editor, model, index):
        value, ok = editor.text().toDouble()
        if editor.validator().state == QValidator.Acceptable:
            model.setData(index, QVariant(value), Qt.DisplayRole)


#-------------------------------------------------------------------------------
# StandarItemModel class to display Coals in a QTableView
#-------------------------------------------------------------------------------

class StandardItemModelCoal(QStandardItemModel):

    def __init__(self, case):
        QStandardItemModel.__init__(self)
        self.headers = [
            self.tr("Coal name"), self.tr("Coal value"), self.tr("Coal unit"),
            self.tr("Coal Temp. \nname"), self.tr("Coal Temp. \nvalue"), self.tr("Coal Temp. \nunit")]
        self.setColumnCount(len(self.headers))
        self.dataCoal = []
        self.__case = case


    def setBoundaryFromLabel(self, label):
        self.modelBoundary = Boundary('coal_inlet', label, self.__case)


    def data(self, index, role):
        if not index.isValid():
            return QVariant()
        if role == Qt.DisplayRole:
            return QVariant(self.dataCoal[index.row()][index.column()])
        return QVariant()


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled
        elif index.column() in [1,4]:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable
        else:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable


    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return QVariant(self.headers[section])
        return QVariant()


    def setData(self, index, value, role):
        row = index.row()
        col = index.column()
        if not hasattr(self, "modelBoundary"):
            log.debug("ERROR in setData (StandardItemModelCoal) : no Boundary model defined")
            return
        v, ok = value.toDouble()
        self.dataCoal[row][col] = v
        if col == 1:
            self.modelBoundary.setCoalFlow(v, row)
        elif col == 4:
            self.modelBoundary.setCoalTemperature(v, row)
        self.emit(SIGNAL("dataChanged(const QModelIndex &, const QModelIndex &)"), index, index)
        return True


    def insertItem(self, nameCoal, valCoal, unitCoal, nameCoalTemp, valCoalTemp, unitCoalTemp):
        line = [nameCoal, valCoal, unitCoal, nameCoalTemp, valCoalTemp, unitCoalTemp]
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
        self.__case = case
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
        self.modelBoundary = Boundary('coal_inlet', label, self.__case)


    def data(self, index, role):
        if not index.isValid():
            return QVariant()
        if role == Qt.DisplayRole:
            classe = index.row()
            coal   = index.column()
            if classe < self.coalClassesNumber[coal]:
                try:
                    return QVariant(self.ratio[coal][classe])
                except:
                    log.debug("ERROR no data for self.ratio[%i][%i] "%(coal, classe))
        return QVariant()


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled
        elif index.row() >= self.coalClassesNumber[index.column()]:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable
        else:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable

    
    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return QVariant("Coal" + " " + str(section+1))
        if orientation == Qt.Vertical and role == Qt.DisplayRole:
            return QVariant("Class" + " " + str(section+1))
        return QVariant()

    
    def setData(self, index, value, role):
        if not hasattr(self, "modelBoundary"):
            log.debug("ERROR in setData (StandardItemModelCoalMass) : no Boundary model defined")
            return
        classe = index.row()
        coal   = index.column()
        v, ok = value.toDouble()
        self.ratio[coal][classe] = v
        log.debug("setData v = %f "%v)

        liste = self.modelBoundary.getCoalRatios(coal)
        lastValue = 0
        for iclasse in range(0, self.coalClassesNumber[coal]-1):
            lastValue += self.ratio[coal][iclasse]

        if lastValue < 100.+ 1e-6 :
            liste[classe] = self.ratio[coal][classe]
            lastValue = 100 - lastValue
            self.ratio[coal][self.coalClassesNumber[coal]-1] = lastValue
            liste[self.coalClassesNumber[coal]-1] = lastValue
            self.modelBoundary.setCoalRatios(coal, liste)
##             self.__getRatioLastClass(coal)
        else :
##             self.ratio[coal][classe].set(model.getClassCoalRatio(coal, classe))
            self.ratio[coal][classe] = liste[classe]

        self.emit(SIGNAL("dataChanged(const QModelIndex &, const QModelIndex &)"), index, index)
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
        
        
    def setup(self, case):
        """
        Setup the widget
        """
        self.setupUi(self)
        self.__case = case
        self.__boundary = None
        
        # Coal table
        self.__modelCoal = StandardItemModelCoal(self.__case)
        self.tableViewCoal.setModel(self.__modelCoal)
        delegateValue = ValueDelegate(self.tableViewCoal)
        self.tableViewCoal.setItemDelegateForColumn(1, delegateValue)
        self.tableViewCoal.setItemDelegateForColumn(4, delegateValue)

        import Pages.CoalCombustionModel as CoalCombustion
        coalModel =  CoalCombustion.CoalCombustionModel(self.__case)
        if coalModel.getCoalCombustionModel() != "off" :
            import Pages.CoalThermoChemistry as CoalThermoChemistry
            coalModel = CoalThermoChemistry.CoalThermoChemistryModel("dp_FCP", self.__case)
            self.__coalNumber = coalModel.getCoals().getNumber()
            self.__coalClassesNumber = []
            for coal in range(0, self.__coalNumber):
                self.__coalClassesNumber.append(coalModel.getCoals().getCoal(coal+1).getClassesNumber())
        else :
            self.__coalNumber = 0
            self.__coalClassesNumber = [0]

        self.__ratio = self.__coalNumber*[0]
        for i in range(0, self.__coalNumber) :
            self.__ratio[i] = self.__coalClassesNumber[i]*[0]

        # Coal mass table
        self.__modelCoalMass = StandardItemModelCoalMass(self.__case, self.__coalNumber, self.__coalClassesNumber)
        self.tableViewCoalMass.setModel(self.__modelCoalMass)

        delegateValueMass = ValueDelegate(self.tableViewCoalMass)
        for c in range(self.__modelCoalMass.columnCount()):
            self.tableViewCoalMass.setItemDelegateForColumn(c, delegateValueMass)

        # Combo models
        self.__modelTypeInlet = ComboModel(self.comboBoxTypeInlet, 2, 1)
        self.__modelTypeInlet.addItem(self.tr("Air"), 'airflow')
        self.__modelTypeInlet.addItem(self.tr("Air & Coal"), 'coalflow')

        self.__modelAirVelocity = ComboModel(self.comboBoxAirVelocity, 6, 1)
        self.__modelAirVelocity.addItem(self.tr("Velocity"), 'norm')
        self.__modelAirVelocity.addItem(self.tr("Mass flow rate"), 'flow1')
        self.__modelAirVelocity.addItem(self.tr("Volumic flow rate"), 'flow2')
        self.__modelAirVelocity.addItem(self.tr("Velocity and direction"), 'norm+direction')
        self.__modelAirVelocity.addItem(self.tr("Mass flow rate and direction"), 'flow1+direction')
        self.__modelAirVelocity.addItem(self.tr("Volumic flow rate and direction"), 'flow2+direction')

        # Validators
        validatorAirFlow = DoubleValidator(self.lineEditAirVelocity)
        validatorTemp = DoubleValidator(self.lineEditTemperature, min=0.)
        validatorXCoal = DoubleValidator(self.lineEditXVelocityCoal)
        validatorYCoal = DoubleValidator(self.lineEditYVelocityCoal)
        validatorZCoal = DoubleValidator(self.lineEditZVelocityCoal)
        
        # Apply validators
        self.lineEditAirVelocity.setValidator(validatorAirFlow)
        self.lineEditTemperature.setValidator(validatorTemp)
        self.lineEditXVelocityCoal.setValidator(validatorXCoal)
        self.lineEditYVelocityCoal.setValidator(validatorYCoal)
        self.lineEditZVelocityCoal.setValidator(validatorZCoal)

# Coals
        self.connect(self.comboBoxTypeInlet,   SIGNAL("activated(const QString&)"),\
                     self.__slotCoalFlowType)
        self.connect(self.comboBoxAirVelocity, SIGNAL("activated(const QString&)"),\
                     self.__slotCoalChoiceVelocity)
        self.connect(self.lineEditAirVelocity, SIGNAL("textChanged(const QString &)"),\
                     self.__slotCoalAirFlow)
        self.connect(self.lineEditTemperature, SIGNAL("textChanged(const QString &)"),\
                     self.__slotTemperature)

#        self.connect(self.comboBoxDirectionVelocityCoal, SIGNAL("activated(const QString&)"), self.__slotCoalDirY)
        self.connect(self.lineEditXVelocityCoal, SIGNAL("textChanged(const QString &)"),\
                     self.__slotCoalDirX)
        self.connect(self.lineEditYVelocityCoal, SIGNAL("textChanged(const QString &)"),\
                     self.__slotCoalDirY)
        self.connect(self.lineEditZVelocityCoal, SIGNAL("textChanged(const QString &)"),\
                     self.__slotCoalDirZ)


    @pyqtSignature("const QString&")
    def __slotCoalFlowType(self, text):
        """
        INPUT inlet type : 'air' or 'air + coal'
        """
        value = self.__modelTypeInlet.dicoV2M[str(text)]
        log.debug("slotCoalFlowType value = %s " % value)
        model = self.__boundary
        label = model.getLabel()
        coal_model = Boundary('coal_inlet', label, self.__case)
        velocity = model.getVelocity()
        self.lineEditAirVelocity.setText(QString(str(velocity)))
        temperature = coal_model.getAirTemperature()
        self.lineEditTemperature.setText(QString(str(temperature)))

        if value == 'airflow':
            self.groupBoxCoal.hide()
            self.groupBoxCoalMass.hide()
            #self.__slotCoalChoiceVelocity()
##                 self.groupBoxDirection.show()
##                 self.groupBoxCoal.show()
        else :
            self.groupBoxCoal.show()
            self.groupBoxCoalMass.show()
            self.groupBoxDirection.show()
            self.groupBoxCoal.show()

            self.__modelCoal.deleteAll()
            self.__modelCoal.setBoundaryFromLabel(label)

            self.__modelCoalMass.deleteAll()
            self.__modelCoalMass.setBoundaryFromLabel(label)

            for coal in range(0, self.__coalNumber):
                # Flow and temperature
                self.__modelCoal.insertItem(self.tr("Coal Flow") + " " + str(coal+1),
                                        coal_model.getCoalFlow(coal), "kg/s",
                                        self.tr("Coal temperature") + " " + str(coal+1),
                                        coal_model.getCoalTemperature(coal), "K")

            for coal in range(0, self.__coalNumber) :
                lastValue = 0.
                for coalClass in range(0, self.__coalClassesNumber[coal]-1):
##                    lastValue += float(coal_model._getClassCoalRatio(coal, coalClass))
                    list = coal_model.getCoalRatios(coal)
                    lastValue += list[coalClass]

                    self.__ratio[coal][coalClass] = list[coalClass]

##                             if (coalClass == self.__coalClassesNumber[coal]-1) :
##                                 self.coale5[coal][coalClass].bind("<<Event>>",TkPage.Callback(self.__getRatioLastClass, coal))
##                             else:
##                                 self.coale5[coal][coalClass].bind("<<Event>>",TkPage.Callback(self.getCoale5, coal, coalClass))
##                                 self.coale5[coal][coalClass].config(fg='black')

                # last class is computed
                coalClass = self.__coalClassesNumber[coal]-1
                lastValue = 100 - lastValue
                self.__ratio[coal][coalClass] = lastValue
                #
                self.__getRatioLastClass(coal)

        self.__modelCoalMass.setRatio(self.__ratio)
        coal_model.setCoalType(value)


    def __getRatioLastClass(self, coal):
        label = self.__boundary.getLabel()
        model = Boundary('coal_inlet', label, self.__case)
##             model.setClassCoalRatio(self.__ratio[coal][self.__coalClassesNumber[coal]-1].get(), coal, self.__coalClassesNumber[coal]-1)



    @pyqtSignature("const QString&")
    def __slotCoalChoiceVelocity(self, text):
        """
        INPUT choice of method of calculation of the velocity for air (coal)
        """
        model = self.__boundary
        label  = model.getLabel()
        choice = model.getVelocityChoice()
        coal_model = Boundary('coal_inlet', label, self.__case)
        type = coal_model.getCoalType()
        new_type = type
        self.type_coal_flow = str(self.comboBoxTypeInlet.currentText())
        new_type = self.type_coal_flow

        coalchoiceflow = self.__modelAirVelocity.dicoV2M[str(self.comboBoxAirVelocity.currentText())]
        log.debug("slotCoalChoiceVelocity coalchoiceflow = %s "%coalchoiceflow)
        if coalchoiceflow != choice:
            new_choice = coalchoiceflow
        else:
            new_choice = choice
        model.setVelocityChoice(new_choice)

 #       self.forgetCoalWindows()
#        self.groupBoxCoal.show()
#        self.__turbulence.groupBoxTurbulence.show()

        self.setWindowsForCoalVelocityChoice(new_choice)
        val = model.getVelocity()
        self.lineEditAirVelocity.setText(QString(str(val)))


    @pyqtSignature("const QString&")
    def __slotCoalAirFlow(self, text):
        self.flow, ok = text.toDouble()
        self.__boundary.setVelocity(self.flow)


    @pyqtSignature("const QString&")
    def __slotTemperature(self, text):
        """
        INPUT air temperature
        """
        label = self.__boundary.getLabel()
        model = Boundary('coal_inlet', label, self.__case)
        temperature, ok = text.toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
            model.setAirTemperature(temperature)


    @pyqtSignature("const QString&")
    def slotCoalDirFlow(self, text):
        """
        INPUT Flow for the velocity
        """
        dico = {"normal" : 0, "vector" : 1, "formula" : 2}
        self.frameVelocity.hide()
        direction = self.modelDirVelocityCoal.dicoV2M[str(text)]
        log.debug("__slotCoalDirY direction = %s "%direction)
        dir = dico[direction]

        model = self.__boundary
        if dir == 1:
            model.updateVelocityChoiceForDirection(dir)
            self.frameVelocity.show()
            x = model.getDirection('direction_x')
            y = model.getDirection('direction_y')
            z = model.getDirection('direction_z')
            self.lineEditXVelocityCoal.setText(QString(str(x)))
            self.lineEditYVelocityCoal.setText(QString(str(y)))
            self.lineEditZVelocityCoal.setText(QString(str(z)))
        else:
            model.updateVelocityChoiceForDirection(dir)
            model.deleteDirectionNodes()


    @pyqtSignature("const QString&")
    def __slotCoalDirX(self, text):
        """
        INPUT value into direction of inlet flow
        """
        model = self.__boundary
        if model.getVelocityChoice()[-9:] == 'direction':
            value, ok = text.toDouble()
            if self.sender().validator().state == QValidator.Acceptable:
                model.setDirection('direction_x', value)
##         else:
##             msg = self.tr("You must select one wall or inlet in the list.")
##             self.stbar.showMessage(msg, 2000)


    @pyqtSignature("const QString&")
    def __slotCoalDirY(self, text):
        """
        INPUT value into direction of inlet flow
        """
        model = self.__boundary
        if model.getVelocityChoice()[-9:] == 'direction':
            value, ok = text.toDouble()
            if self.sender().validator().state == QValidator.Acceptable:
                model.setDirection('direction_y', value)
##         else:
##             msg = self.tr("You must select one wall or inlet in the list.")
##             self.stbar.showMessage(msg, 2000)


    @pyqtSignature("const QString&")
    def __slotCoalDirZ(self, text):
        """
        INPUT value into direction of inlet flow
        """
        model = self.__boundary
        if model.getVelocityChoice()[-9:] == 'direction':
            value, ok = text.toDouble()
            if self.sender().validator().state == QValidator.Acceptable:
                model.setDirection('direction_z', value)
##         else:
##             msg = self.tr("You must select one wall or inlet in the list.")
##             self.stbar.showMessage(msg, 2000)


    def setWindowsForCoalVelocityChoice(self, choice):
        """
        Put windows beyond choice of velocity for inlet nature
        """
        if choice  == 'norm':
            self.labelUnitVelocityAir.setText(QString(str('m/s')))
        elif choice == 'flow1':
            self.labelUnitVelocityAir.setText(QString(str('kg/s')))
        elif choice == 'flow2':
            self.labelUnitVelocityAir.setText(QString(str('m<sup>3</sup>/s')))
    

    def hideWidget(self):
        """
        Hide the widget
        """
        self.groupBoxFlowTemp.hide()
        self.groupBoxCoal.hide()
        self.groupBoxCoalMass.hide()
        self.hide()


    def showWidget(self, b):
        """
        Show the widget
        """
        if self.__coalNumber == 0:
            choice = b.getVelocityChoice()
            label = b.getLabel()
            self.__boundary = Boundary('coal_inlet', label, self.__case)
       
            self.show()
            self.groupBoxFlowTemp.show()
            self.groupBoxCoal.show()
            self.groupBoxCoalMass.show()
    
            type = boundary.getCoalType()
            if type in ("coalflow", "airflow"):
                self.__modelTypeInlet.setItem(str_model=type)
            else:
                msg = "Error :invalid velocity_pressure choice for coal combustion"
                raise ValueError, msg
    
            choice = string.split(choice, '+')[0]
            log.debug("slotSelectBoundary COAL INLET choice = %s "%choice)
            self.__slotCoalFlowType(self.comboBoxTypeInlet.currentText())
    
            self.setWindowsForCoalVelocityChoice(choice)
            self.updateSpecDirWindowForCoal(self.__boundary.getLabel())
        else:
            self.hideWidget()


    def getCoalNumber(self):
        """
        Return the coal number
        """
        return self.__coalNumber


    def updateSpecDirWindowForCoal(self, label):
        """
        Put special window for direction of inlet for coal case
        """
        # ??? Same function as above ???
        #self.updateSpecDirWindow(label)
        model = Boundary('inlet', label, self.__case)
        log.debug("updateSpecDirWindowForCoal model.getVelocityChoice() = %s"%model.getVelocityChoice())
        if model.getVelocityChoice()[-9:] == "direction":
            x = model.getDirection('direction_x')
            y = model.getDirection('direction_y')
            z = model.getDirection('direction_z')
            self.lineEditXVelocityCoal.setText(QString(str(x)))
            self.lineEditYVelocityCoal.setText(QString(str(y)))
            self.lineEditZVelocityCoal.setText(QString(str(z)))
            self.frameCoalVelocity.show()
        else:
            self.frameCoalVelocity.hide()
            
            
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
