# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2018 EDF S.A.
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

from code_saturne.Base.QtCore    import *
from code_saturne.Base.QtGui     import *
from code_saturne.Base.QtWidgets import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------


from code_saturne.Base.Toolbox import GuiParam
from code_saturne.Base.QtPage import IntValidator, DoubleValidator, ComboModel
from code_saturne.Base.QtPage import to_qvariant, from_qvariant, to_text_string

from code_saturne.Pages.LagrangianBoundariesForm import Ui_LagrangianBoundariesForm
from code_saturne.Pages.LocalizationModel import LocalizationModel, Zone
from code_saturne.Pages.LagrangianBoundariesModel import LagrangianBoundariesModel
from code_saturne.Pages.LagrangianModel import LagrangianModel
from code_saturne.Pages.LagrangianStatisticsModel import LagrangianStatisticsModel
from code_saturne.Pages.CoalCombustionModel import CoalCombustionModel


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
        validator = IntValidator(editor, min=0) # nb max sets
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
                        self.tr("Number of sets")]
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
                "inlet" : { "inlet" : self.tr("Particles inlet"),
                            "bounce" : self.tr("Particles rebound"),
                            "outlet" : self.tr("Particles outlet") },
                "outlet" : { "outlet" : self.tr("Particles outlet") },
                "free_inlet_outlet" : { "inlet" : self.tr("Particles inlet"),
                                        "outlet" : self.tr("Particles outlet") },
                "imposed_p_outlet" : { "outlet" : self.tr("Particles outlet") },
                "symmetry" : { "part_symmetry" : self.tr("Particles zero-flux"),
                               "bounce" : self.tr("Particles rebound") }
                }
        else:
            self.dicoM2V = {
                "wall" : { "inlet" : self.tr("Particles inlet"),
                           "bounce" : self.tr("Particles rebound"),
                           "deposit1" : self.tr("Deposition and elimination"),
                           "deposit2" : self.tr("Deposition") },
                "inlet" : { "inlet" : self.tr("Particles inlet"),
                            "bounce" : self.tr("Particles rebound"),
                            "outlet" : self.tr("Particles outlet") },
                "outlet" : { "outlet" : self.tr("Particles outlet") },
                "free_inlet_outlet" : { "inlet" : self.tr("Particles inlet"),
                                        "outlet" : self.tr("Particles outlet") },
                "imposed_p_outlet" : { "outlet" : self.tr("Particles outlet") },
                "symmetry" : { "part_symmetry" : self.tr("Particles zero-flux"),
                               "bounce" : self.tr("Particles rebound")}
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
            n_sets = self.model.getNumberOfSetsValue(label)
            line = [label, nature, interaction, n_sets]
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
            n_sets = from_qvariant(value, int)
            self._data[row][col] = n_sets
            label = self._data[row][0]
            nn = self.model.getNumberOfSetsValue(label)
            label = self._data[row][0]
            self.model.setNumberOfSetsValue(label, n_sets)

        self.dataChanged.emit(index, index)
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
        if QT_API == "PYQT4":
            self.tableViewBoundaries.horizontalHeader().setResizeMode(QHeaderView.Stretch)
        elif QT_API == "PYQT5":
            self.tableViewBoundaries.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)

        delegateInteraction = ParticleBoundaryInteractionDelegate(self.tableViewBoundaries)
        delegateSetNumber = ValueDelegate(self.tableViewBoundaries)
        self.tableViewBoundaries.setItemDelegateForColumn(2,delegateInteraction)
        self.tableViewBoundaries.setItemDelegateForColumn(3,delegateSetNumber)

        self.modelIPOIT = ComboModel(self.comboBoxIPOIT,3,1)
        self.modelIPOIT.addItem(self.tr("Mass flow rate"), "rate")
        self.modelIPOIT.addItem(self.tr("Statistical weight set by values"), "prescribed")

        self.modelIJUVW = ComboModel(self.comboBoxIJUVW,4,1)
        self.modelIJUVW.addItem(self.tr("Fluid velocity"), "fluid")
        self.modelIJUVW.addItem(self.tr("Normal direction velocity"), "norm")
        self.modelIJUVW.addItem(self.tr("Velocity given by values"), "components")

        self.modelIJRTP = ComboModel(self.comboBoxIJRTP,2,1)
        self.modelIJRTP.addItem(self.tr("Fluid temperature"), "fluid")
        self.modelIJRTP.addItem(self.tr("Temperature set by values"), "prescribed")

        self.tableViewBoundaries.clicked[QModelIndex].connect(self.slotSelectBoundary)
        self.modelBoundaries.dataChanged.connect(self.dataChanged)
        self.spinBoxICLAS.valueChanged[int].connect(self.slotICLAS)

        self.lineEditIJNBP.textChanged[str].connect(self.slotIJNBP)
        self.lineEditIJFRE.textChanged[str].connect(self.slotIJFRE)
        self.lineEditICLST.textChanged[str].connect(self.slotICLST)
        self.lineEditIDEBT.textChanged[str].connect(self.slotIDEBT)
        self.comboBoxIPOIT.activated[str].connect(self.slotIPOITChoice)
        self.lineEditIPOIT.textChanged[str].connect(self.slotIPOIT)
        self.lineEditIROPT.textChanged[str].connect(self.slotIROPT)
        self.lineEditIRCOLM.textChanged[str].connect(self.slotIRCOLM)

        self.comboBoxIJUVW.activated[str].connect(self.slotIJUVW)
        self.lineEditIUNO.textChanged[str].connect(self.slotIUNO)
        self.lineEditIUPT.textChanged[str].connect(self.slotIUPT)
        self.lineEditIVPT.textChanged[str].connect(self.slotIVPT)
        self.lineEditIWPT.textChanged[str].connect(self.slotIWPT)

        self.comboBoxIJRTP.activated[str].connect(self.slotIJRTP)
        self.lineEditITPT.textChanged[str].connect(self.slotITPT)
        self.lineEditICPT.textChanged[str].connect(self.slotICPT)
        self.lineEditIEPSI.textChanged[str].connect(self.slotIEPSI)

        self.lineEditIDPT.textChanged[str].connect(self.slotIDPT)
        self.lineEditIVDPT.textChanged[str].connect(self.slotIVDPT)

        self.lineEditINUCHL.textChanged[str].connect(self.slotINUCHL)
        self.lineEditIHPT.textChanged[str].connect(self.slotIHPT)

        # Validators
        validatorIJNBP  = IntValidator(self.lineEditIJNBP, min=0)
        validatorIJFRE  = IntValidator(self.lineEditIJFRE, min=0)
        validatorICLST  = IntValidator(self.lineEditICLST, min=0)
        validatorIDEBT  = DoubleValidator(self.lineEditIDEBT, min=0.)
        validatorIPOIT  = DoubleValidator(self.lineEditIPOIT, min=0.)
        validatorIPOIT.setExclusiveMin(True)
        validatorIROPT  = DoubleValidator(self.lineEditIROPT, min=0.)
        validatorIROPT.setExclusiveMin(True)
        validatorIRCOLM = DoubleValidator(self.lineEditIRCOLM, min=0.)

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
        self.lineEditIRCOLM.setValidator(validatorIRCOLM)

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
        self.groupBoxSetNumber.hide()
        self.groupBoxMain.hide()
        self.groupBoxRate.hide()
        self.groupBoxVelocity.hide()
        self.groupBoxTemperature.hide()
        self.groupBoxDiameter.hide()
        self.groupBoxCoal.hide()


    @pyqtSlot("QModelIndex")
    def slotSelectBoundary(self, index):
        """
        """
        self._hideAllWidgets()
        label, nature, interaction, n_sets = self.modelBoundaries.getItem(index.row())
        self.label = label
        if interaction != "inlet":
            return
        self.model.setCurrentBoundaryNode(nature, label)
        if n_sets > 0:
            self.groupBoxSetNumber.show()
            self.spinBoxICLAS.setMinimum(1)
            self.spinBoxICLAS.setMaximum(n_sets)
            self.spinBoxICLAS.setValue(1)
            self.slotICLAS(1)


    def dataChanged(self, topLeft, bottomRight):
        """
        """
        self.slotSelectBoundary(topLeft)


    @pyqtSlot(int)
    def slotICLAS(self, iset):
        """
        Input ICLAS.
        """
        self.iset = iset
        index = self.tableViewBoundaries.currentIndex()
        label, nature, interaction, n_sets = self.modelBoundaries.getItem(index.row())
        if interaction == "inlet":
            self.model.setCurrentSetNode(self.label, iset)

        self.LM = LagrangianModel(self.case)
        part_model = self.LM.getParticlesModel()

        # Main variables
        self.groupBoxMain.show()
        npart = self.model.getNumberOfParticulesInSetValue(self.label, self.iset)
        self.lineEditIJNBP.setText(str(npart))
        freq = self.model.getInjectionFrequencyValue(self.label, self.iset)
        self.lineEditIJFRE.setText(str(freq))

        self.LSM = LagrangianStatisticsModel(self.case)
        if self.LSM.getGroupOfParticlesValue() > 0:
            igroup = self.model.getParticleGroupNumberValue(self.label, self.iset)
            self.lineEditICLST.setText(str(igroup))
            self.labelICLST.show()
            self.lineEditICLST.show()
        else:
            self.labelICLST.hide()
            self.lineEditICLST.hide()

        # Rate / stat. weight
        self.groupBoxRate.show()
        choice = self.model.getStatisticalWeightChoice(self.label, self.iset)
        self.modelIPOIT.setItem(str_model=choice)
        text = self.modelIPOIT.dicoM2V[choice]
        self.slotIPOITChoice(text)

        # Velocity
        self.groupBoxVelocity.show()
        choice = self.model.getVelocityChoice(self.label, self.iset)
        self.modelIJUVW.setItem(str_model=choice)
        text = self.modelIJUVW.dicoM2V[choice]
        self.slotIJUVW(text)

        # Fouling
        colm = self.model.getFoulingIndexValue(self.label, self.iset)
        self.lineEditIRCOLM.setText(str(colm))

        # Temperature
        status = self.LM.getHeating()
        if part_model == "thermal" and status == "on":
            self.groupBoxTemperature.show()
            choice = self.model.getTemperatureChoice(self.label, self.iset)
            self.modelIJRTP.setItem(str_model=choice)
            text = self.modelIJRTP.dicoM2V[choice]
            self.slotIJRTP(text)

            cp = self.model.getSpecificHeatValue(self.label, self.iset)
            self.lineEditICPT.setText(str(cp))
            eps = self.model.getEmissivityValue(self.label, self.iset)
            self.lineEditIEPSI.setText(str(eps))

        # Coals
        if CoalCombustionModel(self.case).getCoalCombustionModel("only") != 'off':
            self.groupBoxCoal.show()
            icoal = self.model.getCoalNumberValue(self.label, self.iset)
            self.lineEditINUCHL.setText(str(icoal))
            temp  = self.model.getCoalTemperatureValue(self.label, self.iset)
            self.lineEditIHPT.setText(str(temp))

        # Diameter
        self.groupBoxDiameter.show()

        diam = self.model.getDiameterValue(self.label, self.iset)
        vdiam = self.model.getDiameterVarianceValue(self.label, self.iset)
        self.lineEditIDPT.setText(str(diam))
        self.lineEditIVDPT.setText(str(vdiam))

        #Coal
        if CoalCombustionModel(self.case).getCoalCombustionModel("only") != 'off':
            self.labelIROPT.hide()
            self.labelUnitIROPT.hide()
            self.lineEditIROPT.hide()
        else:
            self.labelIROPT.show()
            self.labelUnitIROPT.show()
            self.lineEditIROPT.show()
            rho = self.model.getDensityValue(self.label, self.iset)
            self.lineEditIROPT.setText(str(rho))


    @pyqtSlot(str)
    def slotIJNBP(self, text):
        """
        Input IJNBP.
        """
        if self.lineEditIJNBP.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, int)
            self.model.setNumberOfParticulesInSetValue(self.label, self.iset, value)


    @pyqtSlot(str)
    def slotIJFRE(self, text):
        """
        Input IJFRE.
        """
        if self.lineEditIJFRE.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, int)
            self.model.setInjectionFrequencyValue(self.label, self.iset, value)


    @pyqtSlot(str)
    def slotICLST(self, text):
        """
        Input ICLST.
        """
        if self.lineEditICLST.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, int)
            self.model.setParticleGroupNumberValue(self.label, self.iset, value)


    @pyqtSlot(str)
    def slotIDEBT(self, text):
        """
        Input IDEBT.
        """
        if self.lineEditIDEBT.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.model.setMassFlowRateValue(self.label, self.iset, value)


    @pyqtSlot(str)
    def slotIPOITChoice(self, text):
        """
        Input IPOIT.
        """
        choice = self.modelIPOIT.dicoV2M[str(text)]
        self.model.setStatisticalWeightChoice(self.label, self.iset, choice)
        self.frameMassRate.hide()
        self.frameStatisticalWeight.hide()
        if choice == "rate":
            self.frameMassRate.show()
            rate = self.model.getMassFlowRateValue(self.label, self.iset)
            self.lineEditIDEBT.setText(str(rate))
            self.model.setStatisticalWeightValue(self.label, self.iset, 1)
        elif choice == "prescribed":
            self.frameStatisticalWeight.show()
            weight = self.model.getStatisticalWeightValue(self.label, self.iset)
            self.lineEditIPOIT.setText(str(weight))


    @pyqtSlot(str)
    def slotIPOIT(self, text):
        """
        Input IPOIT.
        """
        if self.lineEditIPOIT.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.model.setStatisticalWeightValue(self.label, self.iset, value)


    @pyqtSlot(str)
    def slotIROPT(self, text):
        """
        Input IROPT.
        """
        if self.lineEditIROPT.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.model.setDensityValue(self.label, self.iset, value)

    @pyqtSlot(str)
    def slotIRCOLM(self, text):
        """
        Input IRCOLM.
        """
        if self.lineEditIRCOLM.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.model.setFoulingIndexValue(self.label, self.iset, value)


    @pyqtSlot(str)
    def slotIJUVW(self, text):
        """
        Input IJUVW.
        """
        choice = self.modelIJUVW.dicoV2M[str(text)]
        self.model.setVelocityChoice(self.label, self.iset, choice)
        self.frameVelocityNorm.hide()
        self.frameVelocityValues.hide()
        if choice == "norm":
            self.frameVelocityNorm.show()
            norm = self.model.getVelocityNormValue(self.label, self.iset)
            self.lineEditIUNO.setText(str(norm))
        elif choice == "components":
            self.frameVelocityValues.show()
            vu = self.model.getVelocityDirectionValue(self.label, self.iset, "x")
            vv = self.model.getVelocityDirectionValue(self.label, self.iset, "y")
            vw = self.model.getVelocityDirectionValue(self.label, self.iset, "z")
            self.lineEditIUPT.setText(str(vu))
            self.lineEditIVPT.setText(str(vv))
            self.lineEditIWPT.setText(str(vw))


    @pyqtSlot(str)
    def slotIUNO(self, text):
        """
        Input IUNO.
        """
        if self.lineEditIUNO.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.model.setVelocityNormValue(self.label, self.iset, value)


    @pyqtSlot(str)
    def slotIUPT(self, text):
        """
        Input IUPT.
        """
        if self.lineEditIUPT.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.model.setVelocityDirectionValue(self.label, self.iset, "x", value)


    @pyqtSlot(str)
    def slotIVPT(self, text):
        """
        Input IVPT.
        """
        if self.lineEditIVPT.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.model.setVelocityDirectionValue(self.label, self.iset, "y", value)


    @pyqtSlot(str)
    def slotIWPT(self, text):
        """
        Input IWPT.
        """
        if self.lineEditIWPT.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.model.setVelocityDirectionValue(self.label, self.iset, "z", value)


    @pyqtSlot(str)
    def slotIJRTP(self, text):
        """
        Input IJRTP.
        """
        choice = self.modelIJRTP.dicoV2M[str(text)]
        self.model.setTemperatureChoice(self.label, self.iset, choice)
        if choice == "prescribed":
            self.frameTemperature.show()
            temp = self.model.getTemperatureValue(self.label, self.iset)
            self.lineEditITPT.setText(str(temp))
        else:
            self.frameTemperature.hide()


    @pyqtSlot(str)
    def slotITPT(self, text):
        """
        Input ITPT.
        """
        if self.lineEditITPT.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.model.setTemperatureValue(self.label, self.iset, value)


    @pyqtSlot(str)
    def slotICPT(self, text):
        """
        Input ICPT.
        """
        if self.lineEditICPT.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.model.setSpecificHeatValue(self.label, self.iset, value)


    @pyqtSlot(str)
    def slotIEPSI(self, text):
        """
        Input IEPSI.
        """
        if self.lineEditIEPSI.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.model.setEmissivityValue(self.label, self.iset, value)


    @pyqtSlot(str)
    def slotIDPT(self, text):
        """
        Input IDPT.
        """
        if self.lineEditIDPT.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.model.setDiameterValue(self.label, self.iset, value)


    @pyqtSlot(str)
    def slotIVDPT(self, text):
        """
        Input IVDPT.
        """
        if self.lineEditIVDPT.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.model.setDiameterVarianceValue(self.label, self.iset, value)


    @pyqtSlot(str)
    def slotINUCHL(self, text):
        """
        Input IHPT.
        """
        if self.lineEditINUCHL.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, int)
            self.model.setCoalNumberValue(self.label, self.iset, value)


    @pyqtSlot(str)
    def slotIHPT(self, text):
        """
        Input IHPT.
        """
        if self.lineEditIHPT.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.model.setCoalTemperatureValue(self.label, self.iset, value)


    def tr(self, text):
        """
        Translation
        """
        return text

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
