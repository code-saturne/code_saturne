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
This module contains the following classes and function:
- MatisseTypeView
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
from Pages.MatisseTypeForm import Ui_MatisseTypeForm
from Pages.MatisseTypeModel import MatisseTypeModel
import Base.QtPage as QtPage

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("MatisseTypeView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class MatisseTypeView(QWidget, Ui_MatisseTypeForm):
    """
    """
    def __init__(self, parent, case):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_MatisseTypeForm.__init__(self)
        self.setupUi(self)

        self.case = case

        self.model = MatisseTypeModel(self.case)

        # Create the Page layout.
        # Combo model
        self.modelType = QtPage.ComboModel(self.comboBoxType,3,1)
        self.modelType.addItem(self.tr("Vault"), 'vault')
        self.modelType.addItem(self.tr("Cask storage"), 'emm')
        self.modelType.addItem(self.tr("Double jacketed wells"), 'djw')

        self.connect(self.comboBoxType,
                     SIGNAL("activated(const QString&)"),
                     self.getMatisseType)

        # initialize

        t = self.model.getMatisseType()
        self.modelType.setItem(str_model=t)


    def getMatisseType(self):
        """
        Input Storage matisse type.
        """
        t = self.modelType.dicoV2M[str(self.comboBoxType.currentText())]
        self.model.setMatisseType(t)


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
