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
- MatisseNetworkView
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
from MatisseNetworkForm import Ui_MatisseNetworkForm
import Base.QtPage as QtPage
import Pages.MatisseTypeModel as MatisseType
import Pages.MatisseGeomModel as MatisseGeom
from Pages.MatisseNetworkModel import MatisseNetworkModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("MatisseNetworkView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class MatisseNetworkView(QWidget, Ui_MatisseNetworkForm):
    """
    """

    def __init__(self, parent, case):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_MatisseNetworkForm.__init__(self)
        self.setupUi(self)

        self.case = case

        self.model = MatisseNetworkModel(self.case)

        model_geom = MatisseGeom.MatisseGeomModel(self.case)
##         self.lineStep = model_geom.getMatisseGeomDoubleVar('ptrres')
##         self.rowStep = model_geom.getMatisseGeomDoubleVar('plgres')
        self.heightStep = model_geom.getMatisseGeomDoubleVar('epchel')
##         self.lineMax = model_geom.getMatisseGeomDoubleVar('nptran')
##         self.rowMax = model_geom.getMatisseGeomDoubleVar('nplgrs')
##         self.heightMax = model_geom.getMatisseGeomDoubleVar('nchest')

        model_mat_type = MatisseType.MatisseTypeModel(self.case)
        self.alveoStat = model_mat_type.node_alveo['status']


        # Create the Page layout.
        self.widgetLine.initWidget(self.case, "network_line")
        self.widgetRow.initWidget(self.case, "network_row")
        # Connections
        self.connect(self.lineEdit_hreso1,
                     SIGNAL("textChanged(const QString &)"),
                     self.getMatisseNetworkVar_hreso1)

        self.connect(self.lineEdit_hreso2,
                     SIGNAL("textChanged(const QString &)"),
                     self.getMatisseNetworkVar_hreso2)

        self.connect(self.lineEdit_hplen,
                     SIGNAL("textChanged(const QString &)"),
                     self.getMatisseNetworkVar_hplen)

        # Validators
        validator_hreso1 = QtPage.DoubleValidator(
            self.lineEdit_hreso1, "validator_hreso1")

        validator_hreso2 = QtPage.DoubleValidator(
            self.lineEdit_hreso2, "validator_hreso2")

        validator_hplen = QtPage.DoubleValidator(
            self.lineEdit_hplen, "validator_hplen")

        self.lineEdit_hreso1.setValidator(validator_hreso1)
        self.lineEdit_hreso2.setValidator(validator_hreso2)
        self.lineEdit_hplen.setValidator(validator_hplen)

        # initialize

        invStat = 'on'
        if self.alveoStat == 'on' :
            invStat = 'off'

        self.nbcellreso1 = self.model.getMatisseNetworkDoubleVar('nbcellreso1')
        self.nbcellreso2 = self.model.getMatisseNetworkDoubleVar('nbcellreso2')
        self.nbcellplen  = self.model.getMatisseNetworkDoubleVar('nbcellplen')

        self.lineEdit_hreso1.setText(QString(str(self.nbcellreso1)))
        self.lineEdit_hreso2.setText(QString(str(self.nbcellreso2)))
        self.lineEdit_hplen.setText(QString(str(self.nbcellplen)))

        self.hreso1 = self.heightStep * self.nbcellreso1
        self.hreso2 = self.heightStep * self.nbcellreso2
        self.hplen  = self.heightStep * self.nbcellplen

        self.model.setMatisseNetworkVar('hreso1', self.hreso1)
        self.model.setMatisseNetworkVar('hreso2', self.hreso2)
        self.model.setMatisseNetworkVar('hplen',  self.hplen)

        self.lineEdit2_1.setText(QString(str(self.hreso1)))
        self.lineEdit2_2.setText(QString(str(self.hreso2)))
        self.lineEdit2_3.setText(QString(str(self.hplen)))


        if invStat == "off":
            self.lineEdit_hreso1.setDisabled(True)

        if self.alveoStat == "off":
            self.lineEdit_hreso2.setDisabled(True)
            self.lineEdit_hplen.setDisabled(True)


        self.lineEdit2_1.setDisabled(True)
        self.lineEdit2_2.setDisabled(True)
        self.lineEdit2_3.setDisabled(True)


    def getMatisseNetworkVar_hreso1(self):
        """
        Input thermiclic load variable.
        """
        self.nbcellreso1, ok = self.lineEdit_hreso1.text().toFloat()
        self.model.setMatisseNetworkVar('nbcellreso1', self.nbcellreso1)
        val = float(self.nbcellreso1) * self.heightStep
        self.model.setMatisseNetworkVar('hreso1', val)


    def getMatisseNetworkVar_hreso2(self):
        """
        Input thermiclic load variable.
        """
        self.nbcellreso2, ok = self.lineEdit_hreso2.text().toFloat()
        self.model.setMatisseNetworkVar('nbcellreso2', self.nbcellreso2)
        val = float(self.nbcellreso2) * self.heightStep
        self.model.setMatisseNetworkVar('hreso2', val)


    def getMatisseNetworkVar_hplen(self):
        """
        Input thermiclic load variable.
        """
        self.nbcellplen, ok = self.lineEdit_hplen.text().toFloat()
        self.model.setMatisseNetworkVar('nbcellplen', self.nbcellplen)
        val = float(self.nbcellplen) * self.heightStep
        self.model.setMatisseNetworkVar('hplen', val)


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