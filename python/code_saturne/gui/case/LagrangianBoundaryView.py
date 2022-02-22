# -*- coding: utf-8 -*-

# -------------------------------------------------------------------------------

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

# -------------------------------------------------------------------------------

"""
This module contains the following classes:
- ValueDelegate
- StandardItemModelBoundaries
- LagrangianBoundariesView
"""

# -------------------------------------------------------------------------------
# Standard modules
# -------------------------------------------------------------------------------
import logging

# -------------------------------------------------------------------------------
# Third-party modules
# -------------------------------------------------------------------------------

from code_saturne.gui.base.QtCore import *
from code_saturne.gui.base.QtGui import *
from code_saturne.gui.base.QtWidgets import *

# -------------------------------------------------------------------------------
# Application modules import
# -------------------------------------------------------------------------------


from code_saturne.model.Common import GuiParam
from code_saturne.gui.base.QtPage import IntValidator, DoubleValidator, ComboModel
from code_saturne.gui.base.QtPage import from_qvariant, to_text_string

from code_saturne.gui.case.LagrangianBoundaryForm import Ui_LagrangianBoundaryForm
from code_saturne.model.LagrangianBoundariesModel import LagrangianBoundariesModel
from code_saturne.model.LagrangianModel import LagrangianModel
from code_saturne.model.LagrangianStatisticsModel import LagrangianStatisticsModel
from code_saturne.model.CoalCombustionModel import CoalCombustionModel

# -------------------------------------------------------------------------------
# log config
# -------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("LagrangianBoundaryView")
log.setLevel(GuiParam.DEBUG)


# -------------------------------------------------------------------------------
# Main class
# -------------------------------------------------------------------------------

class LagrangianBoundaryView(QWidget, Ui_LagrangianBoundaryForm):
    """
    """

    def __init__(self, parent):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_LagrangianBoundaryForm.__init__(self)
        self.setupUi(self)

    def setup(self, case):

        self.case = case
        self.case.undoStopGlobal()
        self.model = LagrangianBoundariesModel(self.case)
        self.zone = None
        self.dicoM2V = {}
        self.dicoV2M = {}

        self._setConnections()
        self._setValidators()

        self.case.undoStartGlobal()

    def _setValidators(self):
        validatorIJNBP = IntValidator(self.lineEditIJNBP, min=0)
        validatorIJFRE = IntValidator(self.lineEditIJFRE, min=0)
        validatorICLST = IntValidator(self.lineEditICLST, min=0)
        validatorIDEBT = DoubleValidator(self.lineEditIDEBT, min=0.)
        validatorIPOIT = DoubleValidator(self.lineEditIPOIT, min=0.)
        validatorIPOIT.setExclusiveMin(True)
        validatorIROPT = DoubleValidator(self.lineEditIROPT, min=0.)
        validatorIROPT.setExclusiveMin(True)
        validatorIRCOLM = DoubleValidator(self.lineEditIRCOLM, min=0.)
        validatorIUNO = DoubleValidator(self.lineEditIUNO)
        validatorIUPT = DoubleValidator(self.lineEditIUPT)
        validatorIVPT = DoubleValidator(self.lineEditIVPT)
        validatorIWPT = DoubleValidator(self.lineEditIWPT)
        validatorITPT = DoubleValidator(self.lineEditITPT)
        validatorICPT = DoubleValidator(self.lineEditICPT)
        validatorIEPSI = DoubleValidator(self.lineEditIEPSI)
        validatorIDPT = DoubleValidator(self.lineEditIDPT, min=0.)
        validatorIVDPT = DoubleValidator(self.lineEditIVDPT)
        validatorINUCHL = IntValidator(self.lineEditINUCHL, min=0)
        validatorIHPT = DoubleValidator(self.lineEditIHPT)
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

    def _setConnections(self):
        self.comboBoxBoundary.activated[str].connect(self.slotSetParticleBoundary)
        self.lineEditNbSets.editingFinished.connect(self.slotNbSets)
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

    def showWidget(self, zone):
        self.zone = zone
        self.show()
        self._hideSubWidgets()
        self._defineV2MDictionary()
        self._fillComboBoxes()
        self._loadCase()

    def _hideSubWidgets(self):
        self.groupBoxNbSets.hide()
        self.groupBoxSetNumber.hide()
        self.groupBoxMain.hide()
        self.groupBoxRate.hide()
        self.groupBoxVelocity.hide()
        self.groupBoxTemperature.hide()
        self.groupBoxDiameter.hide()
        self.groupBoxCoal.hide()

    def _defineV2MDictionary(self):
        if self.model.getFoulingStatus() == "on":
            self.dicoM2V = {
                "wall": {"inlet": self.tr("Particles inlet"),
                         "bounce": self.tr("Particles rebound"),
                         "deposit1": self.tr("Deposition and elimination"),
                         "deposit2": self.tr("Deposition"),
                         "fouling": self.tr("Fouling")},
                "inlet": {"inlet": self.tr("Particles inlet"),
                          "bounce": self.tr("Particles rebound"),
                          "outlet": self.tr("Particles outlet")},
                "outlet": {"outlet": self.tr("Particles outlet")},
                "free_inlet_outlet": {"inlet": self.tr("Particles inlet"),
                                      "outlet": self.tr("Particles outlet")},
                "imposed_p_outlet": {"outlet": self.tr("Particles outlet")},
                "symmetry": {"part_symmetry": self.tr("Particles symmetry"),
                             "bounce": self.tr("Particles rebound")}
            }[self.zone.getNature()]
        else:
            self.dicoM2V = {
                "wall": {"inlet": self.tr("Particles inlet"),
                         "bounce": self.tr("Particles rebound"),
                         "deposit1": self.tr("Deposition and elimination"),
                         "deposit2": self.tr("Deposition")},
                "inlet": {"inlet": self.tr("Particles inlet"),
                          "bounce": self.tr("Particles rebound"),
                          "outlet": self.tr("Particles outlet")},
                "outlet": {"outlet": self.tr("Particles outlet")},
                "free_inlet_outlet": {"inlet": self.tr("Particles inlet"),
                                      "outlet": self.tr("Particles outlet")},
                "imposed_p_outlet": {"outlet": self.tr("Particles outlet")},
                "symmetry": {"part_symmetry": self.tr("Particles symmetry"),
                             "bounce": self.tr("Particles rebound")}
            }[self.zone.getNature()]
        self.dicoV2M = {}
        for k, v in self.dicoM2V.items():
            self.dicoV2M[v] = k

    def _fillComboBoxes(self):
        self.modelParticleBoundary = ComboModel(self.comboBoxBoundary, 1, 1)
        for key, value in self.dicoM2V.items():
            self.modelParticleBoundary.addItem(value, key)
        self.modelIPOIT = ComboModel(self.comboBoxIPOIT, 2, 1)
        self.modelIPOIT.addItem(self.tr("Mass flow rate"), "rate")
        self.modelIPOIT.addItem(self.tr("Statistical weight set by values"), "prescribed")
        self.modelIJUVW = ComboModel(self.comboBoxIJUVW, 3, 1)
        self.modelIJUVW.addItem(self.tr("Fluid velocity"), "fluid")
        self.modelIJUVW.addItem(self.tr("Normal direction velocity"), "norm")
        self.modelIJUVW.addItem(self.tr("Velocity given by values"), "components")
        self.modelIJRTP = ComboModel(self.comboBoxIJRTP, 2, 1)
        self.modelIJRTP.addItem(self.tr("Fluid temperature"), "fluid")
        self.modelIJRTP.addItem(self.tr("Temperature set by values"), "prescribed")

    def _loadCase(self):
        nature = self.zone.getNature()
        label = self.zone.getLabel()
        self.model.setCurrentBoundaryNode(nature, label)
        interaction = self.model.getBoundaryChoice(nature, label)
        if interaction not in self.dicoM2V.keys():
            print("error: BC '" + label + "' (" + nature + ") type '"
                  + interaction + "' requested\n"
                                  "       but model is not active: set to rebound")
            interaction = 'bounce'
            self.model.setBoundaryChoice(nature, label, interaction)
        self.modelParticleBoundary.setItem(str_model=interaction)
        self.slotSetParticleBoundary(
            self.dicoM2V[interaction])  # this is needed because setItem does not trigger comboBox activation
        if interaction == "inlet":
            nb_sets = self.model.getNumberOfSetsValue()
            self.updateInletDisplay(nb_sets)

    def hideWidget(self):
        self.hide()

    def updateInletDisplay(self, nb_sets):
        self.lineEditNbSets.setText(str(nb_sets))
        if int(nb_sets) > 0:
            self.groupBoxSetNumber.show()
            self.spinBoxICLAS.setMinimum(1)
            self.spinBoxICLAS.setMaximum(nb_sets)
            self.spinBoxICLAS.setValue(1)
            self.slotICLAS(1)
        else:
            self.groupBoxSetNumber.hide()
            self.groupBoxMain.hide()
            self.groupBoxRate.hide()
            self.groupBoxVelocity.hide()
            self.groupBoxTemperature.hide()
            self.groupBoxDiameter.hide()
            self.groupBoxCoal.hide()

    @pyqtSlot()
    def slotNbSets(self):
        nb_sets = from_qvariant(self.lineEditNbSets.text(), to_text_string)
        try:
            nb_sets = int(nb_sets)
        except Exception:
            nb_sets = self.model.getNumberOfSetsValue()
            self.lineEditNbSets.setText(str(nb_sets))
            return
        self.model.setNumberOfSetsValue(nb_sets)
        self.updateInletDisplay(nb_sets)

        return

    @pyqtSlot(str)
    def slotSetParticleBoundary(self, interaction):
        interaction = self.dicoV2M[interaction]
        self.model.setBoundaryChoice(self.zone.getNature(), self.zone.getLabel(), interaction)
        if interaction == "inlet":
            self.groupBoxNbSets.show()
        else:
            self.lineEditNbSets.setText("0")
            self.groupBoxNbSets.hide()

    @pyqtSlot(int)
    def slotICLAS(self, iset):
        """
        Input ICLAS.
        """
        self.iset = iset
        label = self.zone.getLabel()
        interaction = self.dicoV2M[str(self.comboBoxBoundary.currentText())]
        if interaction == "inlet":
            self.model.setCurrentSetNode(iset)

        # Main variables
        self.groupBoxMain.show()
        npart = self.model.getNumberOfParticulesInSetValue(label, self.iset)
        self.lineEditIJNBP.setText(str(npart))
        freq = self.model.getInjectionFrequencyValue(label, self.iset)
        self.lineEditIJFRE.setText(str(freq))

        lagStatisticsModel = LagrangianStatisticsModel(self.case)
        if lagStatisticsModel.getGroupOfParticlesValue() > 0:
            igroup = self.model.getParticleGroupNumberValue(label, self.iset)
            self.lineEditICLST.setText(str(igroup))
            self.labelICLST.show()
            self.lineEditICLST.show()
        else:
            self.labelICLST.hide()
            self.lineEditICLST.hide()

        # Rate / stat. weight
        self.groupBoxRate.show()
        choice = self.model.getStatisticalWeightChoice(label, self.iset)
        self.modelIPOIT.setItem(str_model=choice)
        text = self.modelIPOIT.dicoM2V[choice]
        self.slotIPOITChoice(text)

        # Velocity
        self.groupBoxVelocity.show()
        choice = self.model.getVelocityChoice(label, self.iset)
        self.modelIJUVW.setItem(str_model=choice)
        text = self.modelIJUVW.dicoM2V[choice]
        self.slotIJUVW(text)

        # Fouling
        colm = self.model.getFoulingIndexValue(label, self.iset)
        self.lineEditIRCOLM.setText(str(colm))

        # Temperature
        lagModel = LagrangianModel(self.case)
        part_model = lagModel.getParticlesModel()
        status = lagModel.getHeating()
        if part_model == "thermal" and status == "on":
            self.groupBoxTemperature.show()
            choice = self.model.getTemperatureChoice(label, self.iset)
            self.modelIJRTP.setItem(str_model=choice)
            text = self.modelIJRTP.dicoM2V[choice]
            self.slotIJRTP(text)

            cp = self.model.getSpecificHeatValue(label, self.iset)
            self.lineEditICPT.setText(str(cp))
            eps = self.model.getEmissivityValue(label, self.iset)
            self.lineEditIEPSI.setText(str(eps))

        # Coals
        if CoalCombustionModel(self.case).getCoalCombustionModel("only") != 'off':
            self.groupBoxCoal.show()
            icoal = self.model.getCoalNumberValue(label, self.iset)
            self.lineEditINUCHL.setText(str(icoal))
            temp = self.model.getCoalTemperatureValue(label, self.iset)
            self.lineEditIHPT.setText(str(temp))

        # Diameter
        self.groupBoxDiameter.show()

        diam = self.model.getDiameterValue(label, self.iset)
        vdiam = self.model.getDiameterVarianceValue(label, self.iset)
        self.lineEditIDPT.setText(str(diam))
        self.lineEditIVDPT.setText(str(vdiam))

        # Coal
        if CoalCombustionModel(self.case).getCoalCombustionModel("only") != 'off':
            self.labelIROPT.hide()
            self.labelUnitIROPT.hide()
            self.lineEditIROPT.hide()
        else:
            self.labelIROPT.show()
            self.labelUnitIROPT.show()
            self.lineEditIROPT.show()
            rho = self.model.getDensityValue(label, self.iset)
            self.lineEditIROPT.setText(str(rho))

    @pyqtSlot(str)
    def slotIJNBP(self, text):
        """
        Input IJNBP.
        """
        if self.lineEditIJNBP.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, int)
            self.model.setNumberOfParticulesInSetValue(self.zone.getLabel(), self.iset, value)

    @pyqtSlot(str)
    def slotIJFRE(self, text):
        """
        Input IJFRE.
        """
        if self.lineEditIJFRE.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, int)
            self.model.setInjectionFrequencyValue(self.zone.getLabel(), self.iset, value)

    @pyqtSlot(str)
    def slotICLST(self, text):
        """
        Input ICLST.
        """
        if self.lineEditICLST.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, int)
            self.model.setParticleGroupNumberValue(self.zone.getLabel(), self.iset, value)

    @pyqtSlot(str)
    def slotIDEBT(self, text):
        """
        Input IDEBT.
        """
        if self.lineEditIDEBT.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.model.setMassFlowRateValue(self.zone.getLabel(), self.iset, value)

    @pyqtSlot(str)
    def slotIPOITChoice(self, text):
        """
        Input IPOIT.
        """
        choice = self.modelIPOIT.dicoV2M[str(text)]
        self.model.setStatisticalWeightChoice(self.zone.getLabel(), self.iset, choice)
        self.frameMassRate.hide()
        self.frameStatisticalWeight.hide()
        if choice == "rate":
            self.frameMassRate.show()
            rate = self.model.getMassFlowRateValue(self.zone.getLabel(), self.iset)
            self.lineEditIDEBT.setText(str(rate))
            self.model.setStatisticalWeightValue(self.zone.getLabel(), self.iset, 1)
        elif choice == "prescribed":
            self.frameStatisticalWeight.show()
            weight = self.model.getStatisticalWeightValue(self.zone.getLabel(), self.iset)
            self.lineEditIPOIT.setText(str(weight))

    @pyqtSlot(str)
    def slotIPOIT(self, text):
        """
        Input IPOIT.
        """
        if self.lineEditIPOIT.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.model.setStatisticalWeightValue(self.zone.getLabel(), self.iset, value)

    @pyqtSlot(str)
    def slotIROPT(self, text):
        """
        Input IROPT.
        """
        if self.lineEditIROPT.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.model.setDensityValue(self.zone.getLabel(), self.iset, value)

    @pyqtSlot(str)
    def slotIRCOLM(self, text):
        """
        Input IRCOLM.
        """
        if self.lineEditIRCOLM.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.model.setFoulingIndexValue(self.zone.getLabel(), self.iset, value)

    @pyqtSlot(str)
    def slotIJUVW(self, text):
        """
        Input IJUVW.
        """
        choice = self.modelIJUVW.dicoV2M[str(text)]
        self.model.setVelocityChoice(self.zone.getLabel(), self.iset, choice)
        self.frameVelocityNorm.hide()
        self.frameVelocityValues.hide()
        if choice == "norm":
            self.frameVelocityNorm.show()
            norm = self.model.getVelocityNormValue(self.zone.getLabel(), self.iset)
            self.lineEditIUNO.setText(str(norm))
        elif choice == "components":
            self.frameVelocityValues.show()
            vu = self.model.getVelocityDirectionValue(self.zone.getLabel(), self.iset, "x")
            vv = self.model.getVelocityDirectionValue(self.zone.getLabel(), self.iset, "y")
            vw = self.model.getVelocityDirectionValue(self.zone.getLabel(), self.iset, "z")
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
            self.model.setVelocityNormValue(self.zone.getLabel(), self.iset, value)

    @pyqtSlot(str)
    def slotIUPT(self, text):
        """
        Input IUPT.
        """
        if self.lineEditIUPT.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.model.setVelocityDirectionValue(self.zone.getLabel(), self.iset, "x", value)

    @pyqtSlot(str)
    def slotIVPT(self, text):
        """
        Input IVPT.
        """
        if self.lineEditIVPT.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.model.setVelocityDirectionValue(self.zone.getLabel(), self.iset, "y", value)

    @pyqtSlot(str)
    def slotIWPT(self, text):
        """
        Input IWPT.
        """
        if self.lineEditIWPT.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.model.setVelocityDirectionValue(self.zone.getLabel(), self.iset, "z", value)

    @pyqtSlot(str)
    def slotIJRTP(self, text):
        """
        Input IJRTP.
        """
        choice = self.modelIJRTP.dicoV2M[str(text)]
        self.model.setTemperatureChoice(self.zone.getLabel(), self.iset, choice)
        if choice == "prescribed":
            self.frameTemperature.show()
            temp = self.model.getTemperatureValue(self.zone.getLabel(), self.iset)
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
            self.model.setTemperatureValue(self.zone.getLabel(), self.iset, value)

    @pyqtSlot(str)
    def slotICPT(self, text):
        """
        Input ICPT.
        """
        if self.lineEditICPT.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.model.setSpecificHeatValue(self.zone.getLabel(), self.iset, value)

    @pyqtSlot(str)
    def slotIEPSI(self, text):
        """
        Input IEPSI.
        """
        if self.lineEditIEPSI.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.model.setEmissivityValue(self.zone.getLabel(), self.iset, value)

    @pyqtSlot(str)
    def slotIDPT(self, text):
        """
        Input IDPT.
        """
        if self.lineEditIDPT.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.model.setDiameterValue(self.zone.getLabel(), self.iset, value)

    @pyqtSlot(str)
    def slotIVDPT(self, text):
        """
        Input IVDPT.
        """
        if self.lineEditIVDPT.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.model.setDiameterVarianceValue(self.zone.getLabel(), self.iset, value)

    @pyqtSlot(str)
    def slotINUCHL(self, text):
        """
        Input IHPT.
        """
        if self.lineEditINUCHL.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, int)
            self.model.setCoalNumberValue(self.zone.getLabel(), self.iset, value)

    @pyqtSlot(str)
    def slotIHPT(self, text):
        """
        Input IHPT.
        """
        if self.lineEditIHPT.validator().state == QValidator.Acceptable:
            value = from_qvariant(text, float)
            self.model.setCoalTemperatureValue(self.zone.getLabel(), self.iset, value)

# -------------------------------------------------------------------------------
# End
# -------------------------------------------------------------------------------
