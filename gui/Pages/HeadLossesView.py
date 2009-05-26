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
This module defines the HeadLosses model data management.

This module contains the following classes:
- HeadLosses
- HeadLossesView
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

from Base.Toolbox import GuiParam
from HeadLossesForm import Ui_HeadLossesForm
from HeadLossesAdvancedOptionsDialogForm import Ui_HeadLossesAdvancedOptionsDialogForm
import Base.QtPage as QtPage
from Pages.HeadLossesModel import HeadLossesModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("HeadLossesView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Main view class
#-------------------------------------------------------------------------------

class HeadLossesView(QWidget, Ui_HeadLossesForm):
    """
    Class to open HeadLosses Page.
    """
    def __init__(self, parent=None, case=None):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_HeadLossesForm.__init__(self)
        self.setupUi(self)

        self.case = case
        self.model = HeadLossesModel(self.case)

        # Combo model

        self.modelTurbModel = QtPage.ComboModel(self.comboBoxTurbModel,10,1)

        self.modelTurbModel.addItem(self.tr("No model (i.e. laminar flow)"), "off")
        self.modelTurbModel.addItem(self.tr("Mixing length"), "mixing_length")
        self.modelTurbModel.addItem(self.tr("k-epsilon"), "k-epsilon")
        self.modelTurbModel.addItem(self.tr("k-epsilon Linear Production"), "k-epsilon-PL")
        self.modelTurbModel.addItem(self.tr("Rij-epsilon LLR"), "Rij-epsilon")
        self.modelTurbModel.addItem(self.tr("Rij-epsilon SSG"), "Rij-SSG")
        self.modelTurbModel.addItem(self.tr("v2f (phi model)"), "v2f-phi")
        self.modelTurbModel.addItem(self.tr("k-omega SST"), "k-omega-SST")
        self.modelTurbModel.addItem(self.tr("LES (Smagorinsky)"), "LES_Smagorinsky")
        self.modelTurbModel.addItem(self.tr("LES (classical dynamic model)"), "LES_dynamique")

        # Connections

        self.connect(self.comboBoxTurbModel, SIGNAL("activated(const QString&)"), self.slotHeadLossesModel)
        self.connect(self.pushButtonAdvanced, SIGNAL("clicked()"), self.slotAdvancedOptions)
        self.connect(self.lineEditLength, SIGNAL("textChanged(const QString &)"), self.slotLengthScale)

        # Frames display

        self.frameAdvanced.hide()
        self.frameLength.hide()

        # Validator

        validator = QtPage.DoubleValidator(self.lineEditLength, min=0.0)
        validator.setExclusiveMin(True)
        #validator.setFixup(self.model.defaultHeadLossesValues()['length_scale'])
        self.lineEditLength.setValidator(validator)


        # Update the HeadLosses models list with the calculation features

        for turb in self.model.turbModel:
            if turb not in self.model.HeadLossesModelsList():
                self.modelTurbModel.disableItem(str_model=turb)

        # Select the HeadLosses model

        model = self.model.getHeadLossesModel()
        self.modelTurbModel.setItem(str_model=model)
        self.slotHeadLossesModel(self.comboBoxTurbModel.currentText())

        # Length scale

        l_scale = self.model.getLengthScale()
        self.lineEditLength.setText(QString(str(l_scale)))


    @pyqtSignature("const QString&")
    def slotLengthScale(self, text):
        """
        Private slot.
        Input XLOMLG.
        """
        l_scale, ok = text.toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
            self.model.setLengthScale(l_scale)


    @pyqtSignature("const QString&")
    def slotHeadLossesModel(self, text):
        """
        Private slot.
        Input ITURB.
        """
        model = self.modelTurbModel.dicoV2M[str(text)]
        self.model.setHeadLossesModel(model)

        self.frameAdvanced.hide()
        self.frameLength.hide()

        if model == 'mixing_length':
            self.frameLength.show()
            self.frameAdvanced.hide()
            self.model.getLengthScale()
        elif model not in ('off', 'LES_Smagorinsky', 'LES_dynamique'):
            self.frameLength.hide()
            self.frameAdvanced.show()

        if model in ('off', 'LES_Smagorinsky', 'LES_dynamique'):
            self.line.hide()
        else:
            self.line.show()

##         if model in ('LES_Smagorinsky', 'LES_dynamique'):
##             title = self.tr("HeadLosses model")
##             msg   = self.tr("Please report to the informations \n" \
##                             "contain in the user subroutine: 'ussmag'")
##             QMessageBox.warning(self, title, msg)


    @pyqtSignature("")
    def slotAdvancedOptions(self):
        """
        Private slot.
        Ask one popup for advanced specifications
        """
        default = {}
        default['model']         = self.model.getHeadLossesModel()
        default['scale_model']   = self.model.getScaleModel()
        default['gravity_terms'] = self.model.getGravity()
        log.debug("slotAdvancedOptions -> %s" % str(default))

        dialog = HeadLossesAdvancedOptionsDialogView(self, default)
        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotAdvancedOptions -> %s" % str(result))
            self.model.setHeadLossesModel(result['model'])
            self.model.setScaleModel(result['scale_model'])
            self.model.setGravity(result['gravity_terms'])


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