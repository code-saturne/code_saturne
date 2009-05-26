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
- MatisseThermicView
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
from MatisseThermicForm import Ui_MatisseThermicForm
import Base.QtPage as QtPage
from Pages.MatisseThermicModel import MatisseThermicModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("MatisseThermicView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Model class
#-------------------------------------------------------------------------------

class StandardItemModelThermic(QStandardItemModel):

    def __init__(self, case):
        """
        """
        QStandardItemModel.__init__(self)

        self.case = case
        self.setColumnCount(4)
        self.model = MatisseThermicModel(self.case)
        self._initData()
        

    def _initData(self):

        # String Var
        self.imdcnt = ""

        # Double Var
        self.puicon = 0.
        self.tinit  = 0. 
        self.tcrit  = 0.
        self.emicon = 0.
        self.emimur = 0.
        self.hepcnt = 0.
        self.dhpcnt = 0.

        self.variables = [
            ['imdcnt', self.imdcnt],
            ['puicon', self.puicon],
            ['tinit' , self.tinit],
            ['tcrit' , self.tcrit],
            ['emicon', self.emicon],
            ['emimur', self.emimur],
            ['hepcnt', self.hepcnt],
            ['dhpcnt', self.dhpcnt]]

        self.texts = {}
        self.texts['imdcnt'] = (1, self.tr("Natural convection plume modelling"))
        self.texts['puicon'] = (2, self.tr("Canister heating power"), "W")
        self.texts['tinit']  = (3, self.tr("Inlet air temperature"), "<sup>o</sup>C")
        self.texts['tcrit']  = (4, self.tr("Critical outlet air temperature"), "<sup>o</sup>C")
        self.texts['emicon'] = (5, self.tr("Canister emissivity"))
        self.texts['emimur'] = (6, self.tr("Wall emissivity"))
        self.texts['hepcnt'] = (8, self.tr("Natural convection plume erosion height"), "m")
        self.texts['dhpcnt'] = (9, self.tr("Natural convection plume heat flowrate"), "W")

        self.rows_disabled = []
        
        stat = self.model.getNatConvPanacheStatus()
        self.imdcnt = stat
        
        if stat == "off":
            if not 5 in self.rows_disabled : self.rows_disabled.append(5)
            if not 6 in self.rows_disabled : self.rows_disabled.append(6)

        idx = 0
        for variable in self.variables:
            if variable[0] != 'imdcnt':
                val = self.model.getMatisseThermicDoubleVar(variable[0])
                self.variables[idx][1] = val
            idx += 1    

        self.setRowCount(len(self.variables))


    def data(self, index, role):
        if not index.isValid():
            return QVariant()

        if role == Qt.DisplayRole:
            row = index.row()
            var = self.variables[row][0]
            
            if index.column() == 0:
                num = self.texts[var][0]
                return QVariant(num)

            if index.column() == 1:
                txt = self.texts[var][1]
                return QVariant(txt)

            if index.column() == 2:
                val = self.variables[row][1]
                return QVariant(val)

            if index.column() == 3:
                if len(self.texts[var])>2:
                    unit = self.texts[var][2]
                else:
                    unit = ""
                return QVariant(unit)

        if role == Qt.CheckStateRole:
            if index.row() == 0 and index.column() == 2:
                if self.imdcnt == "on":
                    return QVariant(Qt.Checked)
                else:
                    return QVariant(Qt.Unchecked)

        return QVariant()


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled
        if index.row() in self.rows_disabled:
            return Qt.ItemIsSelectable
        if index.row() == 0 and index.column() == 2:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsUserCheckable
        if index.column() in [0,1,3]:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable
        else:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable


    def setData(self, index, value, role):
        row = index.row()

        if index.column() == 2:

            if index.row() == 1: # imdcnt
                v, ok = value.toInt()
                if v == Qt.Unchecked:
                    self.imdcnt = "off"
                    if not 5 in self.rows_disabled : self.rows_disabled.append(5)
                    if not 6 in self.rows_disabled : self.rows_disabled.append(6)
                else:
                    self.imdcnt = "on"
                    if 5 in self.rows_disabled : self.rows_disabled.remove(5)
                    if 6 in self.rows_disabled : self.rows_disabled.remove(6)
                self.model.setNatConvPanacheStatus(self.imdcnt)
            
            else:
                tag = self.variables[row][0]
                num = self.texts[tag][0]
                var = self.variables[row][1]
                v, ok = value.toDouble()
                var = v 
                self.model.setMatisseThermicVar(tag,var)

        self.emit(SIGNAL("dataChanged(const QModelIndex &, const QModelIndex &)"), index, index)
        return True



#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class MatisseThermicView(QWidget, Ui_MatisseThermicForm):
    """
    """
    
    def __init__(self, parent, case):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_MatisseThermicForm.__init__(self)
        self.setupUi(self)
        
        self.case = case

        # Create the Page layout.

        self.modelThermic = StandardItemModelThermic(self.case)
        self.tableView.setModel(self.modelThermic)
        self.tableView.resizeColumnsToContents()
        
        self.widgetLine.initWidget(self.case, "thermal_line")
        self.widgetRow.initWidget(self.case, "thermal_row")
        self.widgetHeight.initWidget(self.case, "thermal_height")

        
##     def _initModel(self):
##         """
##         Instantiate the matisse type modelling class.
##         """
##         self.model = MatisseThermicModel(self.case)
##         self.model_mat_type = MatisseType.MatisseTypeModel(self.case)
##         model_geom = MatisseGeom.MatisseGeomModel(self.case)
##         self.lineStep = model_geom.getMatisseGeomDoubleVar('ptrres')
##         self.rowStep = model_geom.getMatisseGeomDoubleVar('plgres')
##         self.heightStep = model_geom.getMatisseGeomDoubleVar('epchel')

##         self.lineMax = model_geom.getMatisseGeomDoubleVar('nptran')
##         self.rowMax = model_geom.getMatisseGeomDoubleVar('nplgrs')
##         self.heightMax = model_geom.getMatisseGeomDoubleVar('nchest')
        

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