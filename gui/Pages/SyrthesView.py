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
This module contains the following class:
- SyrthesView
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
from SyrthesForm import Ui_SyrthesForm
from Pages.SolutionDomainModel import SolutionDomainModel
import Base.QtPage as QtPage

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("SyrthesView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class SyrthesView(QWidget, Ui_SyrthesForm):
    """
    """
    
    def __init__(self, parent, case):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_SyrthesForm.__init__(self)
        self.setupUi(self)
        
        self.case = case

        # Create the Page layout.

        # Configure Faces Selection Widget
        self.widgetFaces.setTagName("faces_syrthes")
        self.widgetFaces.setCase(self.case)

        # Connections
        self.connect(self.radioButtonCouplingOn,  SIGNAL("clicked()"), self.slotSyrthesCoupling)
        self.connect(self.radioButtonCouplingOff, SIGNAL("clicked()"), self.slotSyrthesCoupling)

        self.connect(self.radioButton2DOn,  SIGNAL("clicked()"), self.slotSyrthes2dMesh)
        self.connect(self.radioButton2DOff, SIGNAL("clicked()"), self.slotSyrthes2dMesh)

        # Signals for the Faces Selection Widget
        self.connect(self.widgetFaces.pushButtonDelete, SIGNAL("clicked()"), self.delListboxSyrthes) 

        # initialize

        if SolutionDomainModel(self.case).getSyrthesCouplingStatus() == 'on':
            self.radioButtonCouplingOn.setChecked(True)
            self.radioButtonCouplingOff.setChecked(False)
            self.groupBox.show()
            if SolutionDomainModel(self.case).getSyrthes2dMeshStatus() == 'on':
                self.radioButton2DOn.setChecked(True)
                self.radioButton2DOff.setChecked(False)
            else:
                self.radioButton2DOn.setChecked(False)
                self.radioButton2DOff.setChecked(True)
        else:
            self.radioButtonCouplingOn.setChecked(False)
            self.radioButtonCouplingOff.setChecked(True)
            self.groupBox.hide()
            SolutionDomainModel(self.case).setSyrthesCouplingStatus('off')

        #for node in self.node_syrthes.xmlGetNodeList('faces_syrthes'):
        result = {}
        result = SolutionDomainModel(self.case).getSyrthesFaces()
        if result:
            self.widgetFaces.insertItemFromDico(result)
        

    @pyqtSignature("")
    def slotSyrthesCoupling(self):
        """
        Do we have a syrthes coupling ?
        """
        if self.radioButtonCouplingOn.isChecked():
            self.syr_on_off = "on"
        else:
            self.syr_on_off = "off"
        answer = self.syr_on_off
        SolutionDomainModel(self.case).setSyrthesCouplingStatus(answer)

        # Hide/Show of the group box is handled by signals in the Form
        if SolutionDomainModel(self.case).getSyrthes2dMeshStatus() == 'on':
            self.radioButton2DOn.setChecked(True)
            self.radioButton2DOff.setChecked(False)
        else:
            self.radioButton2DOn.setChecked(False)
            self.radioButton2DOff.setChecked(True)


    @pyqtSignature("")
    def slotSyrthes2dMesh(self):
        """
        Is the mesh for Syrthes a 2D mesh ?
        This command is invoked with 2 arguments. The first is the name of the
        button subwidget that has toggled. The second is a boolean value
        indicating whether the button subwidget is selected.
        """
        if self.radioButton2DOn.isChecked():
            self.syr2d_on_off  = "on"
        else:
            self.syr2d_on_off  = "off"
        SolutionDomainModel(self.case).setSyrthes2dMeshStatus(self.syr2d_on_off)
        if self.syr2d_on_off == 'off':
##            SolutionDomainModel(self.case).setSyrthes2dMeshStatus('on')
##        else:
##            SolutionDomainModel(self.case).delSyrthes2dMeshNode()
            self.radioButton2DOn.setChecked(False)
            self.radioButton2DOff.setChecked(True)


    @pyqtSignature("")
    def delListboxSyrthes(self):
        """
        Delete the selection from the listView.
        """
        self.widgetFaces.slotDelItem()
        self.syr2d_on_off = "off"
        self.slotSyrthes2dMesh()
        

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