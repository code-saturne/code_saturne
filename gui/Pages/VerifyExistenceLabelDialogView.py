# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2020 EDF S.A.
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

#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------

from code_saturne.Base.QtCore    import *
from code_saturne.Base.QtGui     import *
from code_saturne.Base.QtWidgets import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.Pages.VerifyExistenceLabelDialogForm import Ui_VerifyExistenceLabelDialogForm
from code_saturne.Base.QtPage import RegExpValidator

#-------------------------------------------------------------------------------
# Advanced dialog
#-------------------------------------------------------------------------------

class VerifyExistenceLabelDialogView(QDialog, Ui_VerifyExistenceLabelDialogForm):
    """
    Advanced dialog
    """
    def __init__(self, parent, default):
        """
        Constructor
        """
        QDialog.__init__(self, parent)

        Ui_VerifyExistenceLabelDialogForm.__init__(self)
        self.setupUi(self)

        self.setWindowTitle(self.tr("New label"))
        self.default = default
        self.default['list'].remove(self.default['label'])
        self.result = self.default.copy()

        v = RegExpValidator(self.lineEdit, self.default['regexp'])
        self.lineEdit.setValidator(v)

        self.lineEdit.textChanged[str].connect(self.slotLabel)


    @pyqtSlot(str)
    def slotLabel(self, text):
        """
        Get label.
        """
        if self.sender().validator().state == QValidator.Acceptable:
            label = str(text)

            if label in self.default['list']:
                self.buttonBox.setDisabled(True)
            else:
                self.buttonBox.setDisabled(False)
                self.result['label'] = label
        else:
            self.buttonBox.setDisabled(True)


    def get_result(self):
        """
        Return the result dictionary.
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
        QDialog.reject(self)


    def tr(self, text):
        """
        Translation
        """
        return text

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
