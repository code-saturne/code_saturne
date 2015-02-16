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
This module defines the DarcyLaw model data management.

This module contains the following classes:
- DarcyLawView
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import sys, logging

#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------

from PyQt4.QtCore import *
from PyQt4.QtGui  import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.Base.Toolbox import GuiParam
from code_saturne.Base.QtPage import ComboModel, setGreenColor, DoubleValidator
from code_saturne.Base.QtPage import to_qvariant, from_qvariant, to_text_string
from code_saturne.Pages.DarcyLawForm import Ui_DarcyLawForm
from code_saturne.Pages.LocalizationModel import LocalizationModel, Zone
from code_saturne.Pages.QMeiEditorView import QMeiEditorView
from code_saturne.Pages.DarcyLawModel import DarcyLawModel
from code_saturne.Pages.DarcyModel import DarcyModel
from code_saturne.Pages.DefineUserScalarsModel import DefineUserScalarsModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("DarcyLawView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# StandarItemModel class to display Head Losses Zones in a QTreeView
#-------------------------------------------------------------------------------


class StandardItemModelDarcyLaw(QStandardItemModel):
    def __init__(self):
        QStandardItemModel.__init__(self)
        self.headers = [self.tr("Label"), self.tr("Zone"),
                        self.tr("Selection criteria")]
        self.setColumnCount(len(self.headers))
        self.dataDarcyLawZones = []


    def data(self, index, role):
        if not index.isValid():
            return to_qvariant()
        if role == Qt.DisplayRole:
            return to_qvariant(self.dataDarcyLawZones[index.row()][index.column()])
        return to_qvariant()

    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled
        else:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable

    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return to_qvariant(self.headers[section])
        return to_qvariant()


    def setData(self, index, value, role):
        self.emit(SIGNAL("dataChanged(const QModelIndex &, const QModelIndex &)"), index, index)
        return True


    def insertItem(self, label, name, local):
        line = [label, name, local]
        self.dataDarcyLawZones.append(line)
        row = self.rowCount()
        self.setRowCount(row+1)


    def getItem(self, row):
        return self.dataDarcyLawZones[row]


#-------------------------------------------------------------------------------
# Main view class
#-------------------------------------------------------------------------------

class DarcyLawView(QWidget, Ui_DarcyLawForm):

    def __init__(self, parent, case):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_DarcyLawForm.__init__(self)
        self.setupUi(self)

        self.case = case
        self.case.undoStopGlobal()

        self.mdl = DarcyLawModel(self.case)

        self.list_scalars = []
        self.m_sca = DefineUserScalarsModel(self.case)
        for s in self.m_sca.getUserScalarNameList():
            self.list_scalars.append((s, self.tr("Additional scalar")))

        # Create the Page layout.

        # Model and QTreeView for Head Losses
        self.modelDarcyLaw = StandardItemModelDarcyLaw()
        self.treeView.setModel(self.modelDarcyLaw)

        # Combo model
        self.modelDarcyLawType = ComboModel(self.comboBoxType, 2, 1)
        self.modelDarcyLawType.addItem(self.tr("User law"), 'user')
        self.modelDarcyLawType.addItem(self.tr("Van Genuchten law"), 'VanGenuchten')

        self.modelNameDiff = ComboModel(self.comboBoxNameDiff,1,1)

        self.modelDiff = ComboModel(self.comboBoxDiff, 2, 1)
        self.modelDiff.addItem(self.tr('constant'), 'constant')
        self.modelDiff.addItem(self.tr('variable'), 'variable')

        # Set up validators

        self.lineEditKs.setValidator(DoubleValidator(self.lineEditKs))
        self.lineEditThetas.setValidator(DoubleValidator(self.lineEditThetas))
        self.lineEditThetar.setValidator(DoubleValidator(self.lineEditThetar))
        self.lineEditN.setValidator(DoubleValidator(self.lineEditN))
        self.lineEditL.setValidator(DoubleValidator(self.lineEditL))
        self.lineEditAlpha.setValidator(DoubleValidator(self.lineEditAlpha))
        self.lineEditLongitudinal.setValidator(DoubleValidator(self.lineEditLongitudinal))
        self.lineEditTransverse.setValidator(DoubleValidator(self.lineEditTransverse))

        self.scalar = ""
        scalar_list = self.m_sca.getUserScalarNameList()
        for s in self.m_sca.getScalarsVarianceList():
            if s in scalar_list: scalar_list.remove(s)

        if scalar_list != []:
            self.scalar = scalar_list[0]
            for scalar in scalar_list:
                self.modelNameDiff.addItem(scalar)

        # Connections
        self.connect(self.treeView,             SIGNAL("clicked(const QModelIndex &)"), self.slotSelectDarcyLawZones)
        self.connect(self.comboBoxType,         SIGNAL("activated(const QString&)"),    self.slotDarcyLaw)
        self.connect(self.lineEditKs,           SIGNAL("textChanged(const QString &)"), self.slotKs)
        self.connect(self.lineEditThetas,       SIGNAL("textChanged(const QString &)"), self.slotThetas)
        self.connect(self.lineEditThetar,       SIGNAL("textChanged(const QString &)"), self.slotThetar)
        self.connect(self.lineEditN,            SIGNAL("textChanged(const QString &)"), self.slotN)
        self.connect(self.lineEditL,            SIGNAL("textChanged(const QString &)"), self.slotL)
        self.connect(self.lineEditAlpha,        SIGNAL("textChanged(const QString &)"), self.slotAlpha)
        self.connect(self.lineEditLongitudinal, SIGNAL("textChanged(const QString &)"), self.slotLongitudinal)
        self.connect(self.lineEditTransverse,   SIGNAL("textChanged(const QString &)"), self.slotTransverse)
        self.connect(self.pushButtonUserLaw,    SIGNAL("clicked()"),                    self.slotFormula)
        self.connect(self.comboBoxNameDiff,     SIGNAL("activated(const QString&)"),    self.slotNameDiff)
        self.connect(self.comboBoxDiff,         SIGNAL("activated(const QString&)"),    self.slotStateDiff)
        self.connect(self.pushButtonDiff,       SIGNAL("clicked()"),                    self.slotFormulaDiff)

        # Initialize Widgets

        self.entriesNumber = -1
        d = self.mdl.getNameAndLocalizationZone()
        liste=[]
        liste=list(d.items())
        t=[]
        for t in liste :
            NamLoc=t[1]
            Lab=t[0 ]
            self.modelDarcyLaw.insertItem(Lab, NamLoc[0],NamLoc[1])

        self.forgetStandardWindows()

        self.case.undoStartGlobal()


    @pyqtSignature("const QModelIndex&")
    def slotSelectDarcyLawZones(self, index):
        label, name, local = self.modelDarcyLaw.getItem(index.row())

        if hasattr(self, "modelScalars"): del self.modelScalars
        log.debug("slotSelectDarcyLawZones label %s " % label )
        self.groupBoxType.show()

        choice = self.mdl.getDarcyLawModel(name)
        self.modelDarcyLawType.setItem(str_model=choice)

        if choice == "user":
            self.groupBoxVanGenuchten.hide()
            self.groupBoxUser.show()
            setGreenColor(self.pushButtonUserLaw, True)
        else:
            self.groupBoxVanGenuchten.show()
            self.groupBoxUser.hide()
            self.initializeVanGenuchten(name)

        if DarcyModel(self.case).getDiffusionType() == 'anisotropic':
            self.groupBoxAnisotropicDiffusion.show()
            value = self.mdl.getDiffusionCoefficient(name, "longitudinal")
            self.lineEditLongitudinal.setText(str(value))
            value = self.mdl.getDiffusionCoefficient(name, "transverse")
            self.lineEditTransverse.setText(str(value))

        if self.scalar == "":
            self.groupBoxDiff.hide()
        else :
            self.groupBoxDiff.show()

            diff_choice =  self.mdl.getScalarDiffusivityChoice(self.scalar, name)
            self.modelDiff.setItem(str_model=diff_choice)
            self.modelNameDiff.setItem(str_model=str(self.scalar))
            if diff_choice  != 'variable':
                self.pushButtonDiff.setEnabled(False)
                setGreenColor(self.pushButtonDiff, False)
            else:
                self.pushButtonDiff.setEnabled(True)
                setGreenColor(self.pushButtonDiff, True)

        # force to variable property
        self.comboBoxDiff.setEnabled(False)

        self.entriesNumber = index.row()


    def initializeVanGenuchten(self, name):
        """
        initialize variables for groupBoxVanGenuchten
        """
        value = self.mdl.getValue(name, "ks")
        self.lineEditKs.setText(str(value))
        value = self.mdl.getValue(name, "thetas")
        self.lineEditThetas.setText(str(value))
        value = self.mdl.getValue(name, "thetar")
        self.lineEditThetar.setText(str(value))
        value = self.mdl.getValue(name, "n")
        self.lineEditN.setText(str(value))
        value = self.mdl.getValue(name, "l")
        self.lineEditL.setText(str(value))
        value = self.mdl.getValue(name, "alpha")
        self.lineEditAlpha.setText(str(value))


    def forgetStandardWindows(self):
        """
        For forget standard windows
        """
        self.groupBoxType.hide()
        self.groupBoxUser.hide()
        self.groupBoxVanGenuchten.hide()
        self.groupBoxAnisotropicDiffusion.hide()
        self.groupBoxDiff.hide()


    @pyqtSignature("const QString &")
    def slotDarcyLaw(self, text):
        """
        Method to call 'getState' with correct arguements for 'rho'
        """
        label, name, local = self.modelDarcyLaw.getItem(self.entriesNumber)
        choice = self.modelDarcyLawType.dicoV2M[str(text)]

        self.mdl.setDarcyLawModel(name, choice)

        if choice == "user":
            self.groupBoxVanGenuchten.hide()
            self.groupBoxUser.show()
            setGreenColor(self.pushButtonUserLaw, True)
        else:
            self.groupBoxVanGenuchten.show()
            self.groupBoxUser.hide()
            self.initializeVanGenuchten(name)


    @pyqtSignature("const QString&")
    def slotKs(self, text):
        """
        """
        label, name, local = self.modelDarcyLaw.getItem(self.entriesNumber)
        if self.sender().validator().state == QValidator.Acceptable:
            val = float(text)
            self.mdl.setValue(name, "ks", val)


    @pyqtSignature("const QString&")
    def slotThetas(self, text):
        """
        """
        label, name, local = self.modelDarcyLaw.getItem(self.entriesNumber)
        if self.sender().validator().state == QValidator.Acceptable:
            val = float(text)
            self.mdl.setValue(name, "thetas", val)


    @pyqtSignature("const QString&")
    def slotThetar(self, text):
        """
        """
        label, name, local = self.modelDarcyLaw.getItem(self.entriesNumber)
        if self.sender().validator().state == QValidator.Acceptable:
            val = float(text)
            self.mdl.setValue(name, "thetar", val)


    @pyqtSignature("const QString&")
    def slotN(self, text):
        """
        """
        label, name, local = self.modelDarcyLaw.getItem(self.entriesNumber)
        if self.sender().validator().state == QValidator.Acceptable:
            val = float(text)
            self.mdl.setValue(name, "n", val)


    @pyqtSignature("const QString&")
    def slotL(self, text):
        """
        """
        label, name, local = self.modelDarcyLaw.getItem(self.entriesNumber)
        if self.sender().validator().state == QValidator.Acceptable:
            val = float(text)
            self.mdl.setValue(name, "l", val)


    @pyqtSignature("const QString&")
    def slotAlpha(self, text):
        """
        """
        label, name, local = self.modelDarcyLaw.getItem(self.entriesNumber)
        if self.sender().validator().state == QValidator.Acceptable:
            val = float(text)
            self.mdl.setValue(name, "alpha", val)


    @pyqtSignature("const QString&")
    def slotLongitudinal(self, text):
        """
        """
        label, name, local = self.modelDarcyLaw.getItem(self.entriesNumber)
        if self.sender().validator().state == QValidator.Acceptable:
            val = float(text)
            self.mdl.setDiffusionCoefficient(name, "longitudinal", val)


    @pyqtSignature("const QString&")
    def slotTransverse(self, text):
        """
        """
        label, name, local = self.modelDarcyLaw.getItem(self.entriesNumber)
        if self.sender().validator().state == QValidator.Acceptable:
            val = float(text)
            self.mdl.setDiffusionCoefficient(name, "transverse", val)


    @pyqtSignature("")
    def slotFormula(self):
        """
        User formula for darcy functions
        """
        label, name, local = self.modelDarcyLaw.getItem(self.entriesNumber)

        exp = self.mdl.getDarcyLawFormula(name)

        if exp == None:
            exp = self.getDefaultDarcyLawFormula(choice)

        if DarcyModel(self.case).getPermeabilityType() == 'anisotropic':
            req = [('capacity',     'Capacity'),
                   ('saturation',   'Saturation'),
                   ('permeability[XX]', 'Permeability'),
                   ('permeability[YY]', 'Permeability'),
                   ('permeability[ZZ]', 'Permeability'),
                   ('permeability[XY]', 'Permeability'),
                   ('permeability[XZ]', 'Permeability'),
                   ('permeability[YZ]', 'Permeability')]
        else:
            req = [('capacity',     'Capacity'),
                   ('saturation',   'Saturation'),
                   ('permeability', 'Permeability')]

        exa = """#example: \n""" + self.mdl.getDefaultDarcyLawFormula()

        sym = [('x', 'cell center coordinate'),
               ('y', 'cell center coordinate'),
               ('z', 'cell center coordinate')]

        dialog = QMeiEditorView(self,
                                check_syntax = self.case['package'].get_check_syntax(),
                                expression = exp,
                                required   = req,
                                symbols    = sym,
                                examples   = exa)
        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormula -> %s" % str(result))
            self.mdl.setDarcyLawFormula(name, result)
            setGreenColor(self.sender(), False)


    @pyqtSignature("const QString &")
    def slotNameDiff(self, text):
        """
        Method to set the variance scalar choosed
        """
        label, name, local = self.modelDarcyLaw.getItem(self.entriesNumber)
        choice = self.modelNameDiff.dicoV2M[str(text)]
        log.debug("slotStateDiff -> %s" % (text))
        self.scalar = str(text)

        self.modelDiff.setItem(str_model=self.mdl.getScalarDiffusivityChoice(self.scalar, name))
        setGreenColor(self.pushButtonDiff, True)


    @pyqtSignature("const QString &")
    def slotStateDiff(self, text):
        """
        Method to set diffusion choice for the coefficient
        """
        label, name, local = self.modelDarcyLaw.getItem(self.entriesNumber)
        choice = self.modelDiff.dicoV2M[str(text)]
        log.debug("slotStateDiff -> %s" % (text))

        if choice != 'variable':
            self.pushButtonDiff.setEnabled(False)
            setGreenColor(self.pushButtonDiff, False)
        else:
            self.pushButtonDiff.setEnabled(True)
            setGreenColor(self.pushButtonDiff, True)

        self.mdl.setScalarDiffusivityChoice(self.scalar, name, choice)


    @pyqtSignature("")
    def slotFormulaDiff(self):
        """
        User formula for the diffusion coefficient
        """
        label, namesca, local = self.modelDarcyLaw.getItem(self.entriesNumber)
        name = self.m_sca.getScalarDiffusivityName(self.scalar)
        exp = self.mdl.getDiffFormula(self.scalar, namesca)
        delay_name = str(self.scalar) + "_delay"
        req = [(str(name), str(self.scalar) + ' diffusion coefficient'),
               (delay_name, str(self.scalar)+ ' delay')]
        exa = ''
        sym = [('x', 'cell center coordinate'),
               ('y', 'cell center coordinate'),
               ('z', 'cell center coordinate'),
               ('saturation', 'saturation')]
        sym.append((str(self.scalar),str(self.scalar)))
        dialog = QMeiEditorView(self,
                                check_syntax = self.case['package'].get_check_syntax(),
                                expression = exp,
                                required   = req,
                                symbols    = sym,
                                examples   = exa)
        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormulaDiff -> %s" % str(result))
            self.mdl.setDiffFormula(self.scalar, namesca, result)
            setGreenColor(self.pushButtonDiff, False)


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
