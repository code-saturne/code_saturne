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
This module defines the main application classes for the Qt GUI.
This GUI provides a simple way to display independante pages, in order to put
informations in the XML document, which reflets the treated case.

This module defines the following classes:
- IdView
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import sys

from code_saturne.Base.QtCore    import *
from code_saturne.Base.QtGui     import *
from code_saturne.Base.QtWidgets import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.Base.IdForm import Ui_IdForm

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

class IdView(QWidget, Ui_IdForm):
    """
    Class for the identity dock widget
    """
    def __init__(self):
        """
        Constructor
        """
        QWidget.__init__(self)
        Ui_IdForm.__init__(self)
        self.setupUi(self)


    def setStudyName(self, s):
        """
        Set the study name in the identity dock widget
        """
        self.lineEdit.setText(str(s))

    def setCaseName(self, s):
        """
        Set the case name in the identity dock widget
        """
        self.lineEdit_2.setText(str(s))

    def setXMLFileName(self, s):
        """
        Set the XML file name in the identity dock widget
        """
        self.lineEdit_3.setText(str(s))

    def set(self, study=None, case=None, filename=None):
        """
        Set names in the identity dock widget
        """
        if study is not None:
            self.lineEdit.setText(str(study))

        if case is not None:
           self.lineEdit_2.setText(str(case))

        if filename is not None:
            self.lineEdit_3.setText(str(filename))


if __name__ == "__main__":

    app = QApplication(sys.argv)

    IdView = IdView()
    IdView.setStudyName("toto")
    IdView.setCaseName("tata")
    IdView.setXMLFileName("titi.xml")
    IdView.show()

    sys.exit(app.exec_())
