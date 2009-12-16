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
This module defines the 'Additional user's scalars' page.

This module contains the following classes:
- LabelDelegate
- InitialValueDelegate
- VarianceDelegate
- MinimumDelegate
- MaximumDelegate
- StandardItemModelScalars
- DefineUserScalarsView
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

from DefineUserScalarsForm import Ui_DefineUserScalarsForm

from LocalizationModel import LocalizationModel
from DefineUserScalarsModel import DefineUserScalarsModel

from Base.Common import LABEL_LENGTH_MAX
from Base.Toolbox import GuiParam
from Base.QtPage import ComboModel, DoubleValidator, RegExpValidator

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("DefineUserScalarsView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Line edit delegate for the label
#-------------------------------------------------------------------------------

class LabelDelegate(QItemDelegate):
    """
    Use of a QLineEdit in the table.
    """
    def __init__(self, parent=None):
        QItemDelegate.__init__(self, parent)
        self.parent = parent
        self.old_plabel = ""


    def createEditor(self, parent, option, index):
        editor = QLineEdit(parent)
        self.old_label = ""
        #editor.installEventFilter(self)
        rx = "[_a-zA-Z][_A-Za-z0-9]{1," + str(LABEL_LENGTH_MAX-1) + "}"
        self.regExp = QRegExp(rx)
        v = RegExpValidator(editor, self.regExp)
        editor.setValidator(v)
        return editor


    def setEditorData(self, editor, index):
        value = index.model().data(index, Qt.DisplayRole).toString()
        self.old_plabel = str(value)
        editor.setText(value)


    def setModelData(self, editor, model, index):
        if not editor.isModified():
            return

        if editor.validator().state == QValidator.Acceptable:
            new_plabel = str(editor.text())

            if new_plabel in model.mdl.getScalarLabelsList():
                default = {}
                default['label']  = self.old_plabel
                default['list']   = model.mdl.getScalarLabelsList()
                default['regexp'] = self.regExp
                log.debug("setModelData -> default = %s" % default)
    
                from VerifyExistenceLabelDialogView import VerifyExistenceLabelDialogView
                dialog = VerifyExistenceLabelDialogView(self.parent, default)
                if dialog.exec_():
                    result = dialog.get_result()
                    new_plabel = result['label']
                    log.debug("setModelData -> result = %s" % result)
                else:
                    new_plabel = self.old_plabel

            model.setData(index, QVariant(QString(new_plabel)), Qt.DisplayRole)

#-------------------------------------------------------------------------------
# Line edit delegate for the initial value
#-------------------------------------------------------------------------------

class InitialValueDelegate(QItemDelegate):
    def __init__(self, parent=None):
        super(InitialValueDelegate, self).__init__(parent)
        self.parent = parent


    def createEditor(self, parent, option, index):
        editor = QLineEdit(parent)
        v = DoubleValidator(editor)
        editor.setValidator(v)
        #editor.installEventFilter(self)
        return editor


    def setEditorData(self, editor, index):
        value = index.model().data(index, Qt.DisplayRole).toString()
        editor.setText(value)


    def setModelData(self, editor, model, index):
        if not editor.isModified():
            return
        if editor.validator().state == QValidator.Acceptable:
            value, ok = editor.text().toDouble()
            for idx in self.parent.selectionModel().selectedIndexes():
                if idx.column() == index.column():
                    [label, old_init, vari, mini, maxi] = model.getData(idx)
                    if model.checkInitMinMax(label, value, mini, maxi):
                        model.setData(idx, QVariant(value), Qt.DisplayRole)

#-------------------------------------------------------------------------------
# Combo box delegate for the variance
#-------------------------------------------------------------------------------

class VarianceDelegate(QItemDelegate):
    """
    Use of a combo box in the table.
    """
    def __init__(self, parent):
        super(VarianceDelegate, self).__init__(parent)
        self.parent   = parent


    def createEditor(self, parent, option, index):
        editor = QComboBox(parent)
        self.modelCombo = ComboModel(editor, 1, 1)
        self.modelCombo.addItem(self.tr("No variance"), "no")
        editor.installEventFilter(self)
        return editor


    def setEditorData(self, editor, index):
        sca = index.model().getData(index)[0]
        if index.model().mdl.getVarianceLabelFromScalarLabel(sca):
            return
        
        l1 = index.model().mdl.getScalarLabelsList()
        for s in index.model().mdl.getScalarsVarianceList():
            if s in l1: l1.remove(s)

        if sca in l1: l1.remove(sca)

        for s in l1:
            self.modelCombo.addItem(s, s)


    def setModelData(self, comboBox, model, index):
        txt = str(comboBox.currentText())
        value = self.modelCombo.dicoV2M[txt]
        model.setData(index, QVariant(value), Qt.DisplayRole)


    def tr(self, text):
        """
        Translation
        """
        return text 

#-------------------------------------------------------------------------------
# Line edit delegate for minimum value
#-------------------------------------------------------------------------------

class MinimumDelegate(QItemDelegate):
    def __init__(self, parent=None):
        super(MinimumDelegate, self).__init__(parent)
        self.parent = parent


    def createEditor(self, parent, option, index):
        editor = QLineEdit(parent)
        v = DoubleValidator(editor)
        editor.setValidator(v)
        #editor.installEventFilter(self)
        return editor


    def setEditorData(self, editor, index):
        value = index.model().data(index, Qt.DisplayRole).toString()
        editor.setText(value)


    def setModelData(self, editor, model, index):
        if not editor.isModified():
            return
        if editor.validator().state == QValidator.Acceptable:
            value, ok = editor.text().toDouble()
            for idx in self.parent.selectionModel().selectedIndexes():
                if idx.column() == index.column():
                    [label, init, vari, old_mini, maxi] = model.getData(idx)
                    if model.checkInitMinMax(label, init, value, maxi):
                        model.setData(idx, QVariant(value), Qt.DisplayRole)

#-------------------------------------------------------------------------------
# Line edit delegate for maximum value
#-------------------------------------------------------------------------------

class MaximumDelegate(QItemDelegate):
    def __init__(self, parent=None):
        super(MaximumDelegate, self).__init__(parent)
        self.parent = parent


    def createEditor(self, parent, option, index):
        editor = QLineEdit(parent)
        v = DoubleValidator(editor)
        editor.setValidator(v)
        #editor.installEventFilter(self)
        return editor


    def setEditorData(self, editor, index):
        value = index.model().data(index, Qt.DisplayRole).toString()
        editor.setText(value)


    def setModelData(self, editor, model, index):
        if not editor.isModified():
            return
        if editor.validator().state == QValidator.Acceptable:
            value, ok = editor.text().toDouble()
            for idx in self.parent.selectionModel().selectedIndexes():
                if idx.column() == index.column():
                    [label, init, vari, mini, old_maxi] = model.getData(idx)
                    if model.checkInitMinMax(label, init, mini, value):
                        model.setData(idx, QVariant(value), Qt.DisplayRole)

#-------------------------------------------------------------------------------
# StandarItemModel class
#-------------------------------------------------------------------------------

class StandardItemModelScalars(QStandardItemModel):
    """
    """
    def __init__(self, parent, mdl, zone):
        """
        """
        QStandardItemModel.__init__(self)

        self.headers = [self.tr("Name"),
                        self.tr("Initial\nvalue"),
                        self.tr("Variance\nof scalar"),
                        self.tr("Minimal\nvalue"),
                        self.tr("Maximal\nvalue")]

        self.setColumnCount(len(self.headers))

        self.toolTipRole = [self.tr("Code_Saturne keyword: NSCAUS"),
                            self.tr("Code_Saturne user subroutine: usinv.F"),
                            self.tr("Code_Saturne keyword: ISCAVR"),
                            self.tr("Code_Saturne keyword: SCAMIN"),
                            self.tr("Code_Saturne keyword: SCAMAX")]

        self._data = []
        self._disable = []
        self.parent = parent
        self.mdl  = mdl
        self._zone = zone


    def data(self, index, role):
        if not index.isValid():
            return QVariant()

        row = index.row()
        col = index.column()

        if role == Qt.ToolTipRole:
            return QVariant(self.toolTipRole[col])
        if role == Qt.DisplayRole:
            return QVariant(self._data[row][col])

        return QVariant()


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled
        if (index.row(), index.column()) in self._disable:
            return Qt.ItemIsSelectable
        return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable


    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return QVariant(self.headers[section])
        return QVariant()


    def setData(self, index, value, role):
        if not index.isValid():
            return Qt.ItemIsEnabled

        # Update the row in the table
        row = index.row()
        col = index.column()

        # Label
        if col == 0:
            old_plabel = self._data[row][col]
            new_plabel = str(value.toString())
            self._data[row][col] = new_plabel
            self.mdl.renameScalarLabel(old_plabel, new_plabel)

        # Numerical values: init values, min and max
        elif col == 1:
            self._data[row][col], ok = value.toDouble()

        elif col in [3, 4]:
            if (row,3) not in self._disable:
                self._data[row][col], ok = value.toDouble()

        # Variance
        elif col == 2:
            variance = str(value.toString())
            self._data[row][col] = variance

            if (row,3) in self._disable: self._disable.remove((row,3))
            if (row,4) in self._disable: self._disable.remove((row,4))
            if variance != "no":
                self._data[row][3] = 0.
                self._data[row][4] = self.mdl.defaultScalarValues()['max_value']
                # Disable values in columns 3 and 4
                if (row,3) not in self._disable: self._disable.append((row,3))
                if (row,4) not in self._disable: self._disable.append((row,4))
            else:
                label = self._data[row][0]
                if self.mdl.getScalarType(label) == 'thermal':
                    if (row,2) not in self._disable: self._disable.append((row,2))
                else:
                    if (row,2) in self._disable: self._disable.remove((row,2))

        # Update the XML doc with the new data
        if col != 0:
            [label, init, vari, mini, maxi] = self._data[row]
            self.mdl.setScalarValues(label, self._zone, init, mini, maxi, vari)

        self.emit(SIGNAL("dataChanged(const QModelIndex &, const QModelIndex &)"), index, index)
        return True


    def getData(self, index):
        row = index.row()
        return self._data[row]


    def newItem(self, existing_label=None):
        """
        Add an item in the table view
        """
        row = self.rowCount()

        label = self.mdl.addUserScalar(self._zone, existing_label)
        min = self.mdl.getScalarMinValue(label)
        max = self.mdl.getScalarMaxValue(label)
        ini = self.mdl.getScalarInitialValue(self._zone, label)
        var = self.mdl.getScalarVariance(label)
        if var in ("", "no variance", "no_variance"):
            var = "no"
        scalar = [label, ini, var, min, max]

        self.setRowCount(row+1)
        self._data.append(scalar)

        if self.mdl.getScalarType(label) == 'thermal':
            if (row,2) not in self._disable: self._disable.append((row,2))
        else:
            if (row,2) in self._disable: self._disable.remove((row,2))


    def getItem(self, row):
        """
        Return the values for an item.
        """
        [label, init, vari, mini, maxi] = self._data[row]
        return label, init, vari, mini, maxi


    def deleteItem(self, row):
        """
        Delete the row in the model.
        """
        log.debug("deleteItem row = %i " % row)

        for tuple in [(row,2), (row,3), (row,4)]:
            if tuple in self._disable: self._disable.remove(tuple)

        del self._data[row]
        row = self.rowCount()
        self.setRowCount(row-1)

        for r in range(self.rowCount()):
            if self.mdl.getScalarType(self._data[r][0]) == 'thermal':
                if (r,2) not in self._disable: self._disable.append((r,2))
            else:
                if (r,2) in self._disable: self._disable.remove((r,2))


    def checkInitMinMax(self, label, init, mini, maxi):
        """
        Verify the coherence between init, mini and maxi
        """
        log.debug("checkInitMinMax")
        OK = 1

        if mini > maxi:
            title = self.tr("Information")
            msg = self.tr("The minimal value is greater than the maximal "\
                          "value. Therefore there will be no clipping for the "\
                          "scalar named:\n\n%1").arg(label)
            QMessageBox.information(self.parent, title, msg)
            return OK

        # if min or max is modified, initial value must be checked for all zones
        for z in LocalizationModel('VolumicZone', self.mdl.case).getZones():
            if z.getNature()['initialization'] == "on":
                name = str(z.getCodeNumber())
                ini = self.mdl.getScalarInitialValue(name, label)
                if name == self._zone:
                    ini = init

                if ini < mini or ini > maxi:
                    title = self.tr("Warning")
                    msg = self.tr("The initial value must be set between the "\
                                  "minimal and the maximal value of the scalar "\
                                  "named: %1 for the volume region: %2").arg(label, z.getLabel())
                    QMessageBox.warning(self.parent, title, msg)
                    OK = 0

        return OK

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class DefineUserScalarsView(QWidget, Ui_DefineUserScalarsForm):
    """
    """
    def __init__(self, parent, case, stbar):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_DefineUserScalarsForm.__init__(self)
        self.setupUi(self)

        self.case = case
        self.mdl = DefineUserScalarsModel(self.case)

        # widget layout

        # 0/ Combo box model for zone
        self.modelZone = ComboModel(self.comboBoxZone, 1, 1)

        # Warning: the Volume region 'all_cells' is mandatory, it is exists always
        zones = LocalizationModel('VolumicZone', self.case).getZones()
        for z in zones:
            if z.getNature()['initialization'] == "on":
                label = z.getLabel()
                name = str(z.getCodeNumber())
                self.modelZone.addItem(self.tr(label), name)
                if label == "all_cells":
                    zone_label = label
                    zone_name  = name

        self.modelZone.setItem(str_model = zone_name)

        # 2/ tableView
        self.modelScalars = StandardItemModelScalars(self, self.mdl, zone_name)
        self.table.horizontalHeader().setResizeMode(QHeaderView.Stretch)
#        self.table.setEditTriggers(QAbstractItemView.DoubleClicked)

        # Delegates
        delegateLabel        = LabelDelegate(self.table)
        delegateInitialValue = InitialValueDelegate(self.table)
        delegateVariance     = VarianceDelegate(self.table)
        delegateMinimum      = MinimumDelegate(self.table)
        delegateMaximum      = MaximumDelegate(self.table)

        self.table.setItemDelegateForColumn(0, delegateLabel)
        self.table.setItemDelegateForColumn(1, delegateInitialValue)
        self.table.setItemDelegateForColumn(2, delegateVariance)
        self.table.setItemDelegateForColumn(3, delegateMinimum)
        self.table.setItemDelegateForColumn(4, delegateMaximum)

        # Connections
        self.connect(self.comboBoxZone,     SIGNAL("activated(const QString&)"), self.slotZone)
        self.connect(self.pushButtonNew,    SIGNAL("clicked()"), self.slotAddScalar)
        self.connect(self.pushButtonDelete, SIGNAL("clicked()"), self.slotDeleteScalar)
        self.connect(self.modelScalars,     SIGNAL("dataChanged(const QModelIndex &, const QModelIndex &)"), self.dataChanged)

        # widget initialization
        self.slotZone(zone_label)


    @pyqtSignature("const QString &")
    def slotZone(self, text):
        """
        Get label of volume zone for initialization
        """
        log.debug("slotZone-> text = %s" % text)
        zone_name = self.modelZone.dicoV2M[str(text)]
        log.debug("slotZone-> name = %s" % zone_name)

        self.table.reset()
        self.modelScalars = StandardItemModelScalars(self, self.mdl, zone_name)
        self.table.setModel(self.modelScalars)

        for label in self.mdl.getScalarLabelsList():
            self.modelScalars.newItem(label)


    @pyqtSignature("")
    def slotAddScalar(self):
        """
        Add a new item in the table when the 'Create' button is pushed.
        """
        self.table.clearSelection()
        self.modelScalars.newItem()


    @pyqtSignature("")
    def slotDeleteScalar(self):
        """
        Just delete the current selected entries from the table and
        of course from the XML file.
        """
        list = []
        for index in self.table.selectionModel().selectedRows():
            row = index.row()
            list.append(row)

        list.sort()
        list.reverse()

        for row in list:
            label, init, vari, mini, maxi = self.modelScalars.getItem(row)
            if self.mdl.getScalarType(label) == 'user':
                self.mdl.deleteScalar(label)
#                self.modelScalars.deleteItem(row)

        self.table.clearSelection()
        txt = str(self.comboBoxZone.currentText())
        self.slotZone(txt)


    @pyqtSignature("const QModelIndex &, const QModelIndex &")
    def dataChanged(self, topLeft, bottomRight):
        for row in range(topLeft.row(), bottomRight.row()+1):
            self.tableView.resizeRowToContents(row)
        for col in range(topLeft.column(), bottomRight.column()+1):
            self.tableView.resizeColumnToContents(col)


    def tr(self, text):
        """
        Translation
        """
        return text

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
