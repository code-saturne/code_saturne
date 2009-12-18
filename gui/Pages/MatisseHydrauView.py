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
- MatisseHydrauView
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
from MatisseHydrauForm import Ui_MatisseHydrauForm
import Base.QtPage as QtPage
import Pages.MatisseTypeModel as MatisseType
from MatisseHydrauModel import MatisseHydrauModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("MatisseHydrauView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Model class
#-------------------------------------------------------------------------------

class StandardItemModelHydrau(QStandardItemModel):

    def __init__(self, case):
        """
        """
        QStandardItemModel.__init__(self)

        self.case = case
        self.model = MatisseHydrauModel(self.case)
        self.model_mat_type = MatisseType.MatisseTypeModel(self.case)

        self.setColumnCount(4)
        self._initData()


    def _initData(self):

        # String Var
        self.icofor = ""
        self.iconlg = ""

        #
        # Double Var
        self.debmas = 0.
        self.pdccha = 0.
        self.pdcfch = 0.
        self.dhchea = 0.
        self.sdchea = 0.
        self.pdcche = 0.
        self.pdccch = 0.
        self.dhches = 0.
        self.sdches = 0.
        self.pdcalg = 0.
        self.pdcatv = 0.
        self.argamt = 0.
        self.pdcslg = 0.
        self.pdcstv = 0.
        self.argavl = 0.
        self.amppdc = 0.
        self.dhalve = 0.
        self.dpvent = 0.

        # association between variables and tag

        self.variables = [
            ['icofor', self.icofor],
            ['iconlg', self.iconlg],
            ['debmas', self.debmas],
            ['pdccha', self.pdccha],
            ['pdcfch', self.pdcfch],
            ['dhchea', self.dhchea],
            ['sdchea', self.sdchea],
            ['pdcche', self.pdcche],
            ['pdccch', self.pdccch],
            ['dhches', self.dhches],
            ['sdches', self.sdches],
            ['pdcalg', self.pdcalg],
            ['pdcatv', self.pdcatv],
            ['argamt', self.argamt],
            ['pdcslg', self.pdcslg],
            ['pdcstv', self.pdcstv],
            ['argavl', self.argavl],
            ['amppdc', self.amppdc],
            ['dhalve', self.dhalve],
            ['dpvent', self.dpvent]]

        self.texts = {}
        self.texts['icofor'] = (1,  self.tr("Forced hydraulic air circulation regime"))
        self.texts['iconlg'] = (2,  self.tr("Canister network in row (staggered arrangement otherwise)"))
        self.texts['debmas'] = (3,  self.tr("Forced circulation air flow"), "Kg/s")
        self.texts['pdccha'] = (4,  self.tr("Inlet chimney diffuser headloss"))
        self.texts['pdcfch'] = (5,  self.tr("Inlet chimney filter headloss"))
        self.texts['dhchea'] = (6,  self.tr("Inlet chimney hydraulic diameter"), "m")
        self.texts['sdchea'] = (7,  self.tr("Inlet chimney flow area"), "m<sup>2</sup>")
        self.texts['pdcche'] = (8,  self.tr("Outlet chimney diffuser headloss"))
        self.texts['pdccch'] = (9,  self.tr("Outlet chimney valve headloss"))
        self.texts['dhches'] = (10, self.tr("Outlet chimney hydraulic diameter"), "m")
        self.texts['sdches'] = (11, self.tr("Outlet chimney flow area"), "m<sup>2</sup>")
        self.texts['pdcalg'] = (12, self.tr("Upstream inlet door longitudinal headloss"))
        self.texts['pdcatv'] = (13, self.tr("Upstream inlet door transversal headloss"))
        self.texts['argamt'] = (14, self.tr("Upstream register incline angle (degree)"), "<sup>o</sup>")
        self.texts['pdcslg'] = (16, self.tr("Downstream outlet door longitudinal headloss"))
        self.texts['pdcstv'] = (17, self.tr("Downstream outlet door transversal headloss"))
        self.texts['argavl'] = (18, self.tr("Downstream register incline angle (degree)"), "<sup>o</sup>")
        self.texts['amppdc'] = (20, self.tr("Network headloss amplifying factor"))
        self.texts['dhalve'] = (23, self.tr("Double jacketed wells hydraulic diameter"), "m")
        self.texts['dpvent'] = (24, self.tr("Inlet/outlet atmospheric pressure difference"), "Pa")

        self.rows_disabled = []

        stat = self.model.getConstrainedConvStatus()
        self.icofor = stat
        if not 2 in self.rows_disabled and stat == "off":
            self.rows_disabled.append(2)

        stat = self.model.getInlineContainerNetworkStatus()
        self.iconlg = stat

        i = 0
        for variable in self.variables :
            if variable[0] not in ['icofor', 'iconlg']:
                val = self.model.getMatisseHydrauDoubleVar(variable[0])
                self.variables[i][1] = val

                if variable[0] == 'dhalve' :
                    stat = self.model_mat_type.node_alveo['status']
                    if not i in self.rows_disabled and stat == "off":
                        self.rows_disabled.append(i)
            i += 1

        self.setRowCount(len(self.variables))


    def data(self, index, role):

        if not index.isValid():
            return QVariant()

        if role == Qt.DisplayRole:
            row = index.row()

            if index.column() == 0:
                var = self.variables[row][0]
                num = self.texts[var][0]
                return QVariant(num)

            if index.column() == 1:
                var = self.variables[row][0]
                txt = self.texts[var][1]
                return QVariant(txt)

            if index.column() == 2:
                var = self.variables[row][0]
                val = self.variables[row][1]
                return QVariant(val)

            if index.column() == 3:
                var = self.variables[row][0]
                if len(self.texts[var])>2:
                    unit = self.texts[var][2]
                else:
                    unit = ""
                return QVariant(unit)

        if role == Qt.CheckStateRole:

            if index.row() == 0 and index.column() == 2:
                if self.icofor == "on":
                    return QVariant(Qt.Checked)
                else:
                    return QVariant(Qt.Unchecked)

            if index.row() == 1 and index.column() == 2:
                if self.iconlg == "on":
                    return QVariant(Qt.Checked)
                else:
                    return QVariant(Qt.Unchecked)

        return QVariant()


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled
        if index.row() in self.rows_disabled:
            return Qt.ItemIsSelectable
        if index.row() in [0,1] and index.column()== 2:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsUserCheckable
        if index.column() in [0,1,3]:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable
        else:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable


    def setData(self, index, value, role):

        if index.column() == 2:

            if index.row() == 0: # icofor
                v, ok = value.toInt()
                if v == Qt.Unchecked:
                    self.icofor = "off"
                    if not 2 in self.rows_disabled:
                        self.rows_disabled.append(2)
                else:
                    self.icofor = "on"
                    if 2 in self.rows_disabled:
                        self.rows_disabled.remove(2)
                self.model.setConstrainedConvStatus(self.icofor)

            elif index.row() == 1: # iconlg
                v, ok = value.toInt()
                if v == Qt.Unchecked:
                    self.iconlg = "off"
                else:
                    self.iconlg = "on"
                self.model.setInlineContainerNetworkStatus(self.iconlg)

            else:
                tag = self.variables[index.row()][0]
                num = self.texts[tag][0]
                var = self.variables[index.row()][1]
                v, ok = value.toDouble()
                var = v
                self.model.setMatisseHydrauVar(num, v) # ???

        self.emit(SIGNAL("dataChanged(const QModelIndex &, const QModelIndex &)"), index, index)
        return True


#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class MatisseHydrauView(QWidget, Ui_MatisseHydrauForm):
    """
    """

    def __init__(self, parent, case):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_MatisseHydrauForm.__init__(self)
        self.setupUi(self)

        self.case = case


        # Create the Page layout.

        self.modelHydrau = StandardItemModelHydrau(self.case)
        self.tableView.setModel(self.modelHydrau)
        self.tableView.setAlternatingRowColors(True)
        self.tableView.resizeColumnsToContents()
        self.tableView.setShowGrid(False)


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