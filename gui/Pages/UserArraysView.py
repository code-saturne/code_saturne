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
This module defines the values of reference.

This module contains the following classes and function:
- UserArraysView
"""

#-------------------------------------------------------------------------------
# Library modules import
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
from UserArraysForm import Ui_UserArraysForm
import Base.QtPage as QtPage
from Pages.UserArraysModel import UserArraysModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("UserArraysView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Main View class
#-------------------------------------------------------------------------------

class UserArraysView(QWidget, Ui_UserArraysForm):
    """
    Class to open User arrays Page.
    """
    def __init__(self, parent, case):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_UserArraysForm.__init__(self)
        self.setupUi(self)

        self.case = case
        self.mdl = UserArraysModel(self.case)

        # Connections

        self.connect(self.lineEditICEL, SIGNAL("textChanged(const QString &)"), self.slotICEL)
        self.connect(self.lineEditRCEL, SIGNAL("textChanged(const QString &)"), self.slotRCEL)
        self.connect(self.lineEditIFAC, SIGNAL("textChanged(const QString &)"), self.slotIFAC)
        self.connect(self.lineEditRFAC, SIGNAL("textChanged(const QString &)"), self.slotRFAC)
        self.connect(self.lineEditIFABOR, SIGNAL("textChanged(const QString &)"), self.slotIFABOR)
        self.connect(self.lineEditRFABOR, SIGNAL("textChanged(const QString &)"), self.slotRFABOR)
        self.connect(self.lineEditIDIMLS, SIGNAL("textChanged(const QString &)"), self.slotIDIMLS)
        self.connect(self.lineEditRDIMLS,   SIGNAL("textChanged(const QString &)"), self.slotRDIMLS)

        # Validators

        validatorICEL = QtPage.IntValidator(self.lineEditICEL, min=0)
        validatorRCEL = QtPage.IntValidator(self.lineEditRCEL, min=0)
        validatorIFAC = QtPage.IntValidator(self.lineEditIFAC, min=0)
        validatorRFAC = QtPage.IntValidator(self.lineEditRFAC, min=0)
        validatorIFABOR = QtPage.IntValidator(self.lineEditIFABOR, min=0)
        validatorRFABOR = QtPage.IntValidator(self.lineEditRFABOR, min=0)
        validatorIDIMLS = QtPage.IntValidator(self.lineEditIDIMLS, min=0)
        validatorRDIMLS = QtPage.IntValidator(self.lineEditRDIMLS, min=0)

        self.lineEditICEL.setValidator(validatorICEL)
        self.lineEditRCEL.setValidator(validatorRCEL)
        self.lineEditIFAC.setValidator(validatorIFAC)
        self.lineEditRFAC.setValidator(validatorRFAC)
        self.lineEditIFABOR.setValidator(validatorIFABOR)
        self.lineEditRFABOR.setValidator(validatorRFABOR)
        self.lineEditIDIMLS.setValidator(validatorIDIMLS)
        self.lineEditRDIMLS.setValidator(validatorRDIMLS)

        # Initialization

        icel   = self.mdl.getIntegerNcelet()
        ifac   = self.mdl.getIntegerNfac()
        ifabor = self.mdl.getIntegerNfabor()
        idimls = self.mdl.getIntegerDimless()

        rcel   = self.mdl.getRealNcelet()
        rfac   = self.mdl.getRealNfac()
        rfabor = self.mdl.getRealNfabor()
        rdimls = self.mdl.getRealDimless()

        self.lineEditICEL.setText(QString(str(icel)))
        self.lineEditIFAC.setText(QString(str(ifac)))
        self.lineEditIFABOR.setText(QString(str(ifabor)))
        self.lineEditIDIMLS.setText(QString(str(idimls)))

        self.lineEditRCEL.setText(QString(str(rcel)))
        self.lineEditRFAC.setText(QString(str(rfac)))
        self.lineEditRFABOR.setText(QString(str(rfabor)))
        self.lineEditRDIMLS.setText(QString(str(rdimls)))


    @pyqtSignature("const QString&")
    def slotICEL(self, text):
        """
        Input number of cells with halo for integer array.
        """
        icel, ok = text.toInt()
        if self.sender().validator().state == QValidator.Acceptable:
            self.mdl.setIntegerNcelet(icel)


    @pyqtSignature("const QString&")
    def slotRCEL(self, text):
        """
        Input number of cells with halo for real array.
        """
        rcel, ok = text.toInt()
        if self.sender().validator().state == QValidator.Acceptable:
            self.mdl.setRealNcelet(rcel)


    @pyqtSignature("const QString&")
    def slotIFAC(self, text):
        """
        Input number of internal faces for integer array.
        """
        ifac, ok = text.toInt()
        if self.sender().validator().state == QValidator.Acceptable:
            self.mdl.setIntegerNfac(ifac)


    @pyqtSignature("const QString&")
    def slotRFAC(self, text):
        """
        Input number of internal faces for real array.
        """
        rfac, ok = text.toInt()
        if self.sender().validator().state == QValidator.Acceptable:
            self.mdl.setRealNfac(rfac)


    @pyqtSignature("const QString&")
    def slotIFABOR(self, text):
        """
        Input number of boundary faces for integer array.
        """
        ifabor, ok = text.toInt()
        if self.sender().validator().state == QValidator.Acceptable:
            self.mdl.setIntegerNfabor(ifabor)


    @pyqtSignature("const QString&")
    def slotRFABOR(self, text):
        """
        Input number of boundary faces for real array.
        """
        rfabor, ok = text.toInt()
        if self.sender().validator().state == QValidator.Acceptable:
            self.mdl.setRealNfabor(rfabor)


    @pyqtSignature("const QString&")
    def slotIDIMLS(self, text):
        """
        Input integer value for integer array.
        """
        idimls, ok = text.toInt()
        if self.sender().validator().state == QValidator.Acceptable:
            self.mdl.setIntegerDimless(idimls)


    @pyqtSignature("const QString&")
    def slotRDIMLS(self, text):
        """
        Input integer value for real array.
        """
        rdimls, ok = text.toInt()
        if self.sender().validator().state == QValidator.Acceptable:
            self.mdl.setRealDimless(rdimls)


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