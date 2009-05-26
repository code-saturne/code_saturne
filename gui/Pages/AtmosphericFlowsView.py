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
This module contains the following classes:
- AtmosphericFlowsView
"""

#-------------------------------------------------------------------------------
# Standard modules
#-------------------------------------------------------------------------------

import logging

#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------
from PyQt4.QtCore import pyqtSignature, SIGNAL
from PyQt4.QtGui  import QWidget

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from Pages.AtmosphericFlowsForm    import Ui_AtmosphericFlowsForm
from Pages.AtmosphericFlowsModel   import AtmosphericFlowsModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("AtmosphericFlowsView")
log.setLevel(logging.DEBUG)
      

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
        
        # Define connection
        self.connect(self.checkBoxMeteoData, SIGNAL("clicked(bool)"), 
                     self.slotCheckBoxMeteoData)
        
        # Initialize the widgets        
        isMeteoDataChecked = model.getMeteoDataStatus() == 'on'
        self.checkBoxMeteoData.setChecked( isMeteoDataChecked )


    @pyqtSignature("bool")
    def slotCheckBoxMeteoData(self, checked):
        """
        Called when checkBox state changed
        """
        status = 'off'
        if checked:
            status = 'on'
        self.__model.setMeteoDataStatus(status)


    def tr(self, text):
        """
        Translation.
        """
        return text

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------