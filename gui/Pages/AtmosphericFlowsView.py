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
- AtmosphericFlowsView
"""

#-------------------------------------------------------------------------------
# Standard modules
#-------------------------------------------------------------------------------

import os, logging

#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------

from PyQt4.QtCore import *
from PyQt4.QtGui  import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.Base.Toolbox import GuiParam
from code_saturne.Pages.AtmosphericFlowsForm   import Ui_AtmosphericFlowsForm
from code_saturne.Pages.AtmosphericFlowsModel  import AtmosphericFlowsModel
from code_saturne.Base.QtPage                  import setGreenColor

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("AtmosphericFlowsView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class AtmosphericFlowsView(QWidget, Ui_AtmosphericFlowsForm):
    """
    Main class
    """
    def __init__(self, parent, case):
        """
        Constructor
        """
        # Setup base classes
        QWidget.__init__(self, parent)
        Ui_AtmosphericFlowsForm.__init__(self)
        self.setupUi(self)

        # create model
        model = AtmosphericFlowsModel(case)
        self.__model = model
        self.case = case

        self.case.undoStopGlobal()

        # Define connection
        self.connect(self.checkBoxMeteoData,
                     SIGNAL("clicked(bool)"),
                     self.__slotCheckBoxMeteoData)
        self.connect(self.pushButtonMeteoData,
                     SIGNAL("pressed()"),
                     self.__slotSearchMeteoData)

        # Initialize the widgets
        isMeteoDataChecked = model.getMeteoDataStatus() == 'on'
        self.checkBoxMeteoData.setChecked(isMeteoDataChecked)
        self.labelMeteoFile.setText(str(self.__model.getMeteoDataFileName()))
        self.labelMeteoData.setEnabled(isMeteoDataChecked)
        self.labelMeteoFile.setEnabled(isMeteoDataChecked)

        self.case.undoStartGlobal()


    @pyqtSignature("bool")
    def __slotCheckBoxMeteoData(self, checked):
        """
        Called when checkBox state changed
        """
        status = 'off'
        if checked:
            status = 'on'

        setGreenColor(self.pushButtonMeteoData, checked)
        self.labelMeteoData.setEnabled(checked)
        self.labelMeteoFile.setEnabled(checked)
        self.__model.setMeteoDataStatus(status)


    @pyqtSignature("")
    def __slotSearchMeteoData(self):
        """
        Select a meteorological file of data
        """
        data = self.case['data_path']
        title = self.tr("Meteorological file of data.")
        filetypes = self.tr("Meteo data (*meteo*);;All Files (*)")
        file = QFileDialog.getOpenFileName(self, title, data, filetypes)
        file = str(file)
        if not file:
            return
        file = os.path.basename(file)
        if file not in os.listdir(data):
            title = self.tr("WARNING")
            msg   = self.tr("This selected file is not in the DATA directory")
            QMessageBox.information(self, title, msg)
        else:
            self.labelMeteoFile.setText(str(file))
            self.__model.setMeteoDataFileName(file)
            setGreenColor(self.pushButtonMeteoData, False)


    def tr(self, text):
        """
        Translation.
        """
        return text

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
