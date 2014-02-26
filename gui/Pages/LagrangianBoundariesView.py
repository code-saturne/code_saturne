# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2014 EDF S.A.
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
- ValueDelegate
- StandardItemModelBoundaries
- LagrangianBoundariesView
"""

#-------------------------------------------------------------------------------
# Standard modules
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


from Base.Toolbox import GuiParam
from Base.QtPage import IntValidator, DoubleValidator, ComboModel
from Base.QtPage import to_qvariant, from_qvariant, to_text_string

from Pages.LagrangianBoundariesForm import Ui_LagrangianBoundariesForm
from Pages.LocalizationModel import LocalizationModel, Zone
from Pages.LagrangianBoundariesModel import LagrangianBoundariesModel
from Pages.LagrangianModel import LagrangianModel
from Pages.LagrangianStatisticsModel import LagrangianStatisticsModel
from Pages.CoalCombustionModel import CoalCombustionModel


#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("LagrangianBoundariesView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Line edit delegate with an integere validator
#-------------------------------------------------------------------------------

class ValueDelegate(QItemDelegate):
    def __init__(self, parent=None):
        super(ValueDelegate, self).__init__(parent)
        self.parent = parent

    def createEditor(self, parent, option, index):
        editor = QLineEdit(parent)
        validator = IntValidator(editor, min=0) # nb max classes
        editor.setValidator(validator)
        return editor

    def setEditorData(self, editor, index):
        value = from_qvariant(index.model().data(index, Qt.DisplayRole), to_text_string)
        editor.setText(value)

    def setModelData(self, editor, model, index):
        if editor.validator().state == QValidator.Acceptable:
            value = from_qvariant(editor.text(), float)
            model.setData(index, to_qvariant(value), Qt.DisplayRole)


#-------------------------------------------------------------------------------
# QComboBox delegate for the particle-boundary interaction
#-------------------------------------------------------------------------------


class ParticleBoundaryInteractionDelegate(QItemDelegate):
    """
    Use of a combo box in the table.
    """
    def __init__(self, parent):
        super(ParticleBoundaryInteractionDelegate, self).__init__(parent)
        self.parent = parent


    def createEditor(self, parent, option, index):
        editor = QComboBox(parent)
        self.combo_mdl = ComboModel(editor,1,1)
        nature = index.model()._data[index.row()][1]
        self.dico = index.model().dicoM2V[nature]
        for k, v in list(self.dico.items()):
            self.combo_mdl.addItem(v, k)
        editor.installEventFilter(self)
        editor.setMinimumWidth(100)
        return editor


    def setEditorData(self, comboBox, index):
        row = index.row()
        col = index.column()
        str_model = index.model()._data[row][col]
        self.combo_mdl.setItem(str_model=str_model)


    def setModelData(self, comboBox, model, index):
        txt   = str(comboBox.currentText())
        value = self.combo_mdl.dicoV2M[txt]
        selectionModel = self.parent.selectionModel()
        for idx in selectionModel.selectedIndexes():
            if idx.column() == index.column():
                model.setData(idx, to_qvariant(value), Qt.DisplayRole)


    def tr(self, text):
        return text


#-------------------------------------------------------------------------------
# StandarItemModel class to display boundaries in a QTreeView
#-------------------------------------------------------------------------------


class StandardItemModelBoundaries(QStandardItemModel):
    def __init__(self, case, model):
        QStandardItemModel.__init__(self)
        self.headers = [self.tr("Label"), self.tr("Nature"),
                        self.tr("Particle-boundary\ninteraction"),
                        self.tr("Number of classes")]
        self.setColumnCount(len(self.headers))
        self.case = case
        self.model = model
        self._data = []

        # Corresponding dict for the nature of the boundary. used in combo delegate.
        if self.model.getFoulingStatus() == "on":
            self.dicoM2V = {
                "wall" : { "inlet" : self.tr("Particles inlet"),
                           "bounce" : self.tr("Particles rebound"),
                           "deposit1" : self.tr("Deposition and elimination"),
                           "deposit2" : self.tr("Deposition"),
                           "fouling" : self.tr("Fouling") },
                "inlet" : { "inlet" : self.tr("Particles inlet") },
                "outlet" : { "outlet" : self.tr("Particles outlet") },
                "symmetry" : { "part_symmetry" : self.tr("Particles zero-flux") }
                }
        else:
            self.dicoM2V = {
                "wall" : { "inlet" : self.tr("Particles inlet"),
                           "bounce" : self.tr("Particles rebound"),
                           "deposit1" : self.tr("Deposition and elimination"),
                           "deposit2" : self.tr("Deposition") },
                "inlet" : { "inlet" : self.tr("Particles inlet") },
                "outlet" : { "outlet" : self.tr("Particles outlet") },
                "symmetry" : { "part_symmetry" : self.tr("Particles zero-flux") }
                }

        self.dicoV2M = {}
        for key in list(self.dicoM2V.keys()):
            dico = self.dicoM2V[key]
            self.dicoV2M[key] = {}
            for k, v in list(dico.items()):
                self.dicoV2M[key][v] = k

        # Initialization
        for zone in LocalizationModel('BoundaryZone', self.case).getZones():
            label = zone.getLabel()
            nature = zone.getNature()
            interaction = self.model.getBoundaryChoice(nature, label)
            nclasses = self.model.getNumberOfClassesValue(label)
            line = [label, nature, interaction, nclasses]
            self._data.append(line)
            row = self.rowCount()
            self.setRowCount(row+1)


    def data(self, index, role):
        if not index.isValid():
            return to_qvariant()

        if role == Qt.DisplayRole:
            row = index.row()
            col = index.column()
            if col == 2:
                nature = self._data[row][1]
                dico = self.dicoM2V[nature]
                return to_qvariant(dico[self._data[row][col]])
            else:
                return to_qvariant(self._data[row][col])

        if role == Qt.ToolTipRole:
            if index.column() == 2:
                return to_qvariant(self.tr("Code_Saturne keyword: IUSCLB"))
            elif index.column() == 3:
                return to_qvariant(self.tr("Code_Saturne keyword: NBCLAS"))
        return to_qvariant()


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled
        elif index.column() in [0,1]:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable
        elif index.column() == 2:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable
        elif index.column() == 3:
            if self._data[index.row()][2] == "inlet":
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

        if col == 2:
            interaction = from_qvariant(value, to_text_string)
            self._data[row][col] = interaction
            label = self._data[row][0]
            nature = self._data[row][1]
            self.model.setBoundaryChoice(nature, label, interaction)
            if nature != "inlet":
                self._data[row][3] = 0

        elif col == 3:
            nclasses = from_qvariant(value, int)
            self._data[row][col] = nclasses
            label = self._data[row][0]
            nn = self.model.getNumberOfClassesValue(label)
            label = self._data[row][0]
            self.model.setNumberOfClassesValue(label, nclasses)

        self.emit(SIGNAL("dataChanged(const QModelIndex &, const QModelIndex &)"), index, index)
        return True


    def getItem(self, row):
        return self._data[row]


    def tr(self, text):
        return text


#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------


class LagrangianBoundariesView(QWidget, Ui_LagrangianBoundariesForm):
    """
    """
    def __init__(self, parent, case):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_LagrangianBoundariesForm.__init__(self)
        self.setupUi(self)

        self.case = case
        self.case.undoStopGlobal()
        self.model = LagrangianBoundariesModel(self.case)

        self.modelBoundaries = StandardItemModelBoundaries(self.case, self.model)
        self.tableViewBoundaries.setModel(self.modelBoundaries)
        self.tableViewBoundaries.setAlternatingRowColors(True)
        self.tableViewBoundaries.horizontalHeader().setResizeMode(QHeaderView.Stretch)

        delegateInteraction = ParticleBoundaryInteractionDelegate(self.tableViewBoundaries)
        delegateClassNumber = ValueDelegate(self.tableViewBoundaries)
        self.tableViewBoundaries.setItemDelegateForColumn(2,delegateInteraction)
        self.tableViewBoundaries.setItemDelegateForColumn(3,delegateClassNumber)

        self.modelIPOIT = ComboModel(self.comboBoxIPOIT,3,1)
        self.modelIPOIT.addItem(self.tr("Volumic flow rate"), "rate")
        self.modelIPOIT.addItem(self.tr("Statistical weight set by values"), "prescribed")
        self.modelIPOIT.addItem(self.tr("User defined statistical weight"), "subroutine")

        self.modelIJUVW = ComboModel(self.comboBoxIJUVW,4,1)
        self.modelIJUVW.addItem(self.tr("Fluid velocity"), "fluid")
        self.modelIJUVW.addItem(self.tr("Normal direction velocity"), "norm")
        self.modelIJUVW.addItem(self.tr("Velocity given by values"), "components")
        self.modelIJUVW.addItem(self.tr("User defined velocity"), "subroutine")

        self.modelIJRTP = ComboModel(self.comboBoxIJRTP,2,1)
        self.modelIJRTP.addItem(self.tr("Temperature set by values"), "prescribed")
        self.modelIJRTP.addItem(self.tr("User defined temperature"), "subroutine")

        self.modelIJRDP = ComboModel(self.comboBoxIJRDP,2,1)
        self.modelIJRDP.addItem(self.tr("Diameter set by values"), "prescribed")
        self.modelIJRDP.addItem(self.tr("User defined diameter"), "subroutine")

        self.modelIRAWCL = ComboModel(self.comboBoxIRAWCL,2,1)
        self.modelIRAWCL.addItem(self.tr("Raw coal"), "raw_coal_as_received")
        self.modelIRAWCL.addItem(self.tr("User defined"), "subroutine")

        self.connect(self.tableViewBoundaries, SIGNAL("clicked(const QModelIndex &)"), self.slotSelectBoundary)
        self.connect(self.modelBoundaries,     SIGNAL("dataChanged(const QModelIndex &, const QModelIndex &)"), self.dataChanged)
        self.connect(self.spinBoxICLAS, SIGNAL("valueChanged(int)"), self.slotICLAS)

        self.connect(self.lineEditIJNBP,  SIGNAL("textChanged(const QString &)"), self.slotIJNBP)
        self.connect(self.lineEditIJFRE,  SIGNAL("textChanged(const QString &)"), self.slotIJFRE)
        self.connect(self.lineEditICLST,  SIGNAL("textChanged(const QString &)"), self.slotICLST)
        self.connect(self.lineEditIDEBT,  SIGNAL("textChanged(const QString &)"), self.slotIDEBT)
        self.connect(self.comboBoxIPOIT,  SIGNAL("activated(const QString&)"), self.slotIPOITChoice)
        self.connect(self.lineEditIPOIT,  SIGNAL("textChanged(const QString &)"), self.slotIPOIT)
        self.connect(self.lineEditIROPT,  SIGNAL("textChanged(const QString &)"), self.slotIROPT)

        self.connect(self.comboBoxIJUVW, SIGNAL("activated(const QString&)"),    self.slotIJUVW)
        self.connect(self.lineEditIUNO,  SIGNAL("textChanged(const QString &)"), self.slotIUNO)
        self.connect(self.lineEditIUPT,  SIGNAL("textChanged(const QString &)"), self.slotIUPT)
        self.connect(self.lineEditIVPT,  SIGNAL("textChanged(const QString &)"), self.slotIVPT)
        self.connect(self.lineEditIWPT,  SIGNAL("textChanged(const QString &)"), self.slotIWPT)

        self.connect(self.comboBoxIJRTP, SIGNAL("activated(const QString&)"),    self.slotIJRTP)
        self.connect(self.lineEditITPT,  SIGNAL("textChanged(const QString &)"), self.slotITPT)
        self.connect(self.lineEditICPT,  SIGNAL("textChanged(const QString &)"), self.slotICPT)
        self.connect(self.lineEditIEPSI, SIGNAL("textChanged(const QString &)"), self.slotIEPSI)

        self.connect(self.comboBoxIJRDP, SIGNAL("activated(const QString&)"),    self.slotIJRDP)
        self.connect(self.lineEditIDPT,  SIGNAL("textChanged(const QString &)"), self.slotIDPT)
        self.connect(self.lineEditIVDPT, SIGNAL("textChanged(const QString &)"), self.slotIVDPT)

        self.connect(self.lineEditINUCHL, SIGNAL("textChanged(const QString &)"), self.slotINUCHL)
        self.connect(self.lineEditIHPT,   SIGNAL("textChanged(const QString &)"), self.slotIHPT)
        self.connect(self.comboBoxIRAWCL, SIGNAL("activated(const QString&)"),    self.slotIRAWCL)

        # Validators
        validatorIJNBP  = IntValidator(self.lineEditIJNBP, min=0)
        validatorIJFRE  = IntValidator(self.lineEditIJFRE, min=0)
        validatorICLST  = IntValidator(self.lineEditICLST, min=0)
        validatorIDEBT  = DoubleValidator(self.lineEditIDEBT, min=0.)
        validatorIPOIT  = DoubleValidator(self.lineEditIPOIT, min=0.)
        validatorIPOIT.setExclusiveMin(True)
        validatorIROPT  = DoubleValidator(self.lineEditIROPT, min=0.)
        validatorIROPT.setExclusiveMin(True)

        validatorIUNO = DoubleValidator(self.lineEditIUNO)
        validatorIUPT = DoubleValidator(self.lineEditIUPT)
        validatorIVPT = DoubleValidator(self.lineEditIVPT)
        validatorIWPT = DoubleValidator(self.lineEditIWPT)

        validatorITPT  = DoubleValidator(self.lineEditITPT)
        validatorICPT  = DoubleValidator(self.lineEditICPT)
        validatorIEPSI = DoubleValidator(self.lineEditIEPSI)

        validatorIDPT  = DoubleValidator(self.lineEditIDPT, min=0.)
        validatorIVDPT = DoubleValidator(self.lineEditIVDPT)

        validatorINUCHL = IntValidator(self.lineEditINUCHL, min=0)
        validatorIHPT   = DoubleValidator(self.lineEditIHPT)

        self.lineEditIJNBP.setValidator(validatorIJNBP)
        self.lineEditIJFRE.setValidator(validatorIJFRE)
        self.lineEditICLST.setValidator(validatorICLST)
        self.lineEditIDEBT.setValidator(validatorIDEBT)
        self.lineEditIPOIT.setValidator(validatorIPOIT)
        self.lineEditIROPT.setValidator(validatorIROPT)

        self.lineEditIUNO.setValidator(validatorIUNO)
        self.lineEditIUPT.setValidator(validatorIUPT)
        self.lineEditIVPT.setValidator(validatorIVPT)
        self.lineEditIWPT.setValidator(validatorIWPT)

        self.lineEditITPT.setValidator(validatorITPT)
        self.lineEditICPT.setValidator(validatorICPT)
        self.lineEditIEPSI.setValidator(validatorIEPSI)

        self.lineEditIDPT.setValidator(validatorIDPT)
        self.lineEditIVDPT.setValidator(validatorIVDPT)

        self.lineEditINUCHL.setValidator(validatorINUCHL)
        self.lineEditIHPT.setValidator(validatorIHPT)

        self._hideAllWidgets()

        self.case.undoStartGlobal()


    def _hideAllWidgets(self):
        self.groupBoxClassNumber.hide()
        self.groupBoxMain.hide()
        self.groupBoxRate.hide()
        self.groupBoxVelocity.hide()
        self.groupBoxTemperature.hide()
        self.groupBoxDiameter.hide()
        self.groupBoxCoal.hide()


    @pyqtSignature("const QModelIndex&")
    def slotSelectBoundary(self, index):
        """
        """
        self._hideAllWidgets()
        label, nature, interaction, nclasses = self.modelBoundaries.getItem(index.row())
        self.label = label
        if interaction != "inlet":
            return
        self.model.setCurrentBoundaryNode(nature, label)
        if nclasses > 0:
            self.groupBoxClassNumber.show()
            self.spinBoxICLAS.setMinimum(1)
            self.spinBoxICLAS.setMaximum(nclasses)
            self.spinBoxICLAS.setValue(1)
            self.slotICLAS(1)


    @pyqtSignature("const QModelIndex &, const QModelIndex &")
    def dataChanged(self, topLeft, bottomRight):
        """
        """
        self.slotSelectBoundary(topLeft)


    @pyqtSignature("int")
    def slotICLAS(self, iclass):
        """
        Input ICLAS.
        """
        self.iclass = iclass
        index = self.tableViewBoundaries.currentIndex()
        label, nature, interaction, nclasses = self.modelBoundaries.getItem(index.row())
        if interaction == "inlet":
            self.model.setCurrentClassNode(self.label, iclass)

        self.LM = LagrangianModel(self.case)
        part_model = self.LM.getParticlesModel()

        # Main variables
        self.groupBoxMain.show()
        npart = self.model.getNumberOfParticulesInClassValue(self.label, self.iclass)
        self.lineEditIJNBP.setText(str(npart))
        freq = self.model.getInjectionFrequencyValue(self.label, self.iclass)
        self.lineEditIJFRE.setText(str(freq))

        self.LSM = LagrangianStatisticsModel(self.case)
        if self.LSM.getGroupOfParticlesValue() > 0:
            igroup = self.model.getParticleGroupNumberValue(self.label, self.iclass)
            self.lineEditICLST.setText(str(igroup))
        else:
            self.labelICLST.setDisabled(True)
            self.lineEditICLST.setDisabled(True)

        # Rate / stat. weight
        self.groupBoxRate.show()
        choice = self.model.getStatisticalWeightChoice(self.label, self.iclass)
        self.modelIPOIT.setItem(str_model=choice)
        text = self.modelIPOIT.dicoM2V[choice]
        self.slotIPOITChoice(from_qvariant(text, to_text_string))

        # Velocity
        self.groupBoxVelocity.show()
        choice = self.model.getVelocityChoice(self.label, self.iclass)
        self.modelIJUVW.setItem(str_model=choice)
        text = self.modelIJUVW.dicoM2V[choice]
        self.slotIJUVW(from_qvariant(text, to_text_string))

        # Temperature
        status = self.LM.getHeating()
        if part_model == "thermal" and status == "on":
            self.groupBoxTemperature.show()
            choice = self.model.getTemperatureChoice(self.label, self.iclass)
            self.modelIJRTP.setItem(str_model=choice)
            text = self.modelIJRTP.dicoM2V[choice]
            self.slotIJRTP(from_qvariant(text, to_text_string))

            cp = self.model.getSpecificHeatValue(self.label, self.iclass)
            self.lineEditICPT.setText(str(cp))
            eps = self.model.getEmissivityValue(self.label, self.iclass)
            self.lineEditIEPSI.setText(str(eps))

        # Coals
        if CoalCombustionModel(self.case).getCoalCombustionModel() != 'off':
            self.groupBoxCoal.show()
            icoal = self.model.getCoalNumberValue(self.label, self.iclass)
            self.lineEditINUCHL.setText(str(icoal))
            temp  = self.model.getCoalTemperatureValue(self.label, self.iclass)
            self.lineEditIHPT.setText(str(temp))
            choice = self.model.getCoalCompositionChoice(self.label, self.iclass)
            self.modelIRAWCL.setItem(str_model=choice)

        # Diameter
        self.groupBoxDiameter.show()
        choice = self.model.getDiameterChoice(self.label, self.iclass)

        if part_model == "coal":
            self.modelIJRDP.setItem(str_model="prescribed")
        else:
            self.modelIJRDP.setItem(str_model=choice)
            text = self.modelIJRDP.dicoM2V[choice]
            self.slotIJRDP(from_qvariant(text, to_text_string))

        if choice == "prescribed":
            self.frameDiameter.show()
            diam = self.model.getDiameterValue(self.label, self.iclass)
            vdiam = self.model.getDiameterVarianceValue(self.label, self.iclass)
            self.lineEditIDPT.setText(str(diam))
            self.lineEditIVDPT.setText(str(vdiam))
        elif choice == "subroutine":
            self.frameDiameter.hide()

        #Coal
        if CoalCombustionModel(self.case).getCoalCombustionModel() != 'off':
            self.labelIROPT.hide()
            self.labelUnitIROPT.hide()
            self.lineEditIROPT.hide()
        else:
            self.labelIROPT.show()
            self.labelUnitIROPT.show()
            self.lineEditIROPT.show()
            rho = self.model.getDensityValue(self.label, self.iclass)
            self.lineEditIROPT.setText(str(rho))


    @pyqtSignature("const QString&")
    def slotIJNBP(self, text):
        """
        Input IJNBP.
        """
        if self.sender().validator().state == QValidator.Acceptable:
            value = from_qvariant(text, int)
            self.model.setNumberOfParticulesInClassValue(self.label, self.iclass, value)


    @pyqtSignature("const QString&")
    def slotIJFRE(self, text):
        """
        Input IJFRE.
        """
        if self.sender().validator().state == QValidator.Acceptable:
            value = from_qvariant(text, int)
            self.model.setInjectionFrequencyValue(self.label, self.iclass, value)


    @pyqtSignature("const QString&")
    def slotICLST(self, text):
        """
        Input ICLST.
        """
        if self.sender().validator().state == QValidator.Acceptable:
            value = from_qvariant(text, int)
            self.model.setParticleGroupNumberValue(self.label, self.iclass, value)


    @pyqtSignature("const QString&")
    def slotIDEBT(self, text):
        """
        Input IDEBT.
        """
        if self.sender().validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.model.setMassFlowRateValue(self.label, self.iclass, value)


    @pyqtSignature("const QString&")
    def slotIPOITChoice(self, text):
        """
        Input IPOIT.
        """
        choice = self.modelIPOIT.dicoV2M[str(text)]
        self.model.setStatisticalWeightChoice(self.label, self.iclass, choice)
        self.frameVolumicRate.hide()
        self.frameStatisticalWeight.hide()
        if choice == "rate":
            self.frameVolumicRate.show()
            rate = self.model.getMassFlowRateValue(self.label, self.iclass)
            self.lineEditIDEBT.setText(str(rate))
            self.model.setStatisticalWeightValue(self.label, self.iclass, 1)
        elif choice == "prescribed":
            self.frameStatisticalWeight.show()
            weight = self.model.getStatisticalWeightValue(self.label, self.iclass)
            self.lineEditIPOIT.setText(str(weight))
        elif choice == "subroutine":
            pass


    @pyqtSignature("const QString&")
    def slotIPOIT(self, text):
        """
        Input IPOIT.
        """
        if self.sender().validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.model.setStatisticalWeightValue(self.label, self.iclass, value)


    @pyqtSignature("const QString&")
    def slotIROPT(self, text):
        """
        Input IROPT.
        """
        if self.sender().validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.model.setDensityValue(self.label, self.iclass, value)


    @pyqtSignature("const QString&")
    def slotIJUVW(self, text):
        """
        Input IJUVW.
        """
        choice = self.modelIJUVW.dicoV2M[str(text)]
        self.model.setVelocityChoice(self.label, self.iclass, choice)
        self.frameVelocityNorm.hide()
        self.frameVelocityValues.hide()
        if choice == "norm":
            self.frameVelocityNorm.show()
            norm = self.model.getVelocityNormValue(self.label, self.iclass)
            self.lineEditIUNO.setText(str(norm))
        elif choice == "components":
            self.frameVelocityValues.show()
            vu = self.model.getVelocityDirectionValue(self.label, self.iclass, "u")
            vv = self.model.getVelocityDirectionValue(self.label, self.iclass, "v")
            vw = self.model.getVelocityDirectionValue(self.label, self.iclass, "w")
            self.lineEditIUPT.setText(str(vu))
            self.lineEditIVPT.setText(str(vv))
            self.lineEditIWPT.setText(str(vw))


    @pyqtSignature("const QString&")
    def slotIUNO(self, text):
        """
        Input IUNO.
        """
        if self.sender().validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.model.setVelocityNormValue(self.label, self.iclass, value)


    @pyqtSignature("const QString&")
    def slotIUPT(self, text):
        """
        Input IUPT.
        """
        if self.sender().validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.model.setVelocityDirectionValue(self.label, self.iclass, "u", value)


    @pyqtSignature("const QString&")
    def slotIVPT(self, text):
        """
        Input IVPT.
        """
        if self.sender().validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.model.setVelocityDirectionValue(self.label, self.iclass, "v", value)


    @pyqtSignature("const QString&")
    def slotIWPT(self, text):
        """
        Input IWPT.
        """
        if self.sender().validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.model.setVelocityDirectionValue(self.label, self.iclass, "w", value)


    @pyqtSignature("const QString&")
    def slotIJRTP(self, text):
        """
        Input IJRTP.
        """
        choice = self.modelIJRTP.dicoV2M[str(text)]
        self.model.setTemperatureChoice(self.label, self.iclass, choice)
        if choice == "prescribed":
            self.frameTemperature.show()
            temp = self.model.getTemperatureValue(self.label, self.iclass)
            self.lineEditITPT.setText(str(temp))
        elif choice == "subroutine":
            self.frameTemperature.hide()


    @pyqtSignature("const QString&")
    def slotITPT(self, text):
        """
        Input ITPT.
        """
        if self.sender().validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.model.setTemperatureValue(self.label, self.iclass, value)


    @pyqtSignature("const QString&")
    def slotICPT(self, text):
        """
        Input ICPT.
        """
        if self.sender().validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.model.setSpecificHeatValue(self.label, self.iclass, value)


    @pyqtSignature("const QString&")
    def slotIEPSI(self, text):
        """
        Input IEPSI.
        """
        if self.sender().validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.model.setEmissivityValue(self.label, self.iclass, value)


    @pyqtSignature("const QString&")
    def slotIJRDP(self, text):
        """
        Input IJRDP.
        """
        choice = self.modelIJRDP.dicoV2M[str(text)]
        self.model.setDiameterChoice(self.label, self.iclass, choice)
        if choice == "prescribed":
            self.frameDiameter.show()
            diam = self.model.getDiameterValue(self.label, self.iclass)
            vdiam = self.model.getDiameterVarianceValue(self.label, self.iclass)
            self.lineEditIDPT.setText(str(diam))
            self.lineEditIVDPT.setText(str(vdiam))
        elif choice == "subroutine":
            self.frameDiameter.hide()


    @pyqtSignature("const QString&")
    def slotIDPT(self, text):
        """
        Input IDPT.
        """
        if self.sender().validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.model.setDiameterValue(self.label, self.iclass, value)


    @pyqtSignature("const QString&")
    def slotIVDPT(self, text):
        """
        Input IVDPT.
        """
        if self.sender().validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.model.setDiameterVarianceValue(self.label, self.iclass, value)


    @pyqtSignature("const QString&")
    def slotINUCHL(self, text):
        """
        Input IHPT.
        """
        if self.sender().validator().state == QValidator.Acceptable:
            value = from_qvariant(text, int)
            self.model.setCoalNumberValue(self.label, self.iclass, value)


    @pyqtSignature("const QString&")
    def slotIHPT(self, text):
        """
        Input IHPT.
        """
        if self.sender().validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.model.setCoalTemperatureValue(self.label, self.iclass, value)


    @pyqtSignature("const QString&")
    def slotIRAWCL(self, text):
        """
        Input IJRDP.
        """
        choice = self.modelIRAWCL.dicoV2M[str(text)]
        self.model.setCoalCompositionChoice(self.label, self.iclass, choice)


    def tr(self, text):
        """
        Translation
        """
        return text

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
