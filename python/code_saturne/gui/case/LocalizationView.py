# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2022 EDF S.A.
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
- VolumicZoneAdvancedView
- LabelDelegate
- CodeNumberDelegate
- BoundaryNatureDelegate
- StandardItemVolumeNature
- FlagBox
- VolumeNatureDelegate
- LocalizationSelectorDelegate
- LocalizationView
- VolumeLocalizationView
- BoundaryLocalizationView
"""

#-------------------------------------------------------------------------------
# Standard modules
#-------------------------------------------------------------------------------

import logging

#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------

from code_saturne.gui.base.QtCore    import *
from code_saturne.gui.base.QtGui     import *
from code_saturne.gui.base.QtWidgets import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.model.Common import LABEL_LENGTH_MAX, GuiParam, GuiLabelManager
from code_saturne.gui.base.QtPage import IntValidator, RegExpValidator
from code_saturne.gui.base.QtPage import from_qvariant, to_text_string
from code_saturne.gui.case.LocalizationForm import Ui_LocalizationForm
from code_saturne.gui.case.VolumicZoneAdvancedDialogForm import Ui_VolumicZoneAdvancedDialogForm
from code_saturne.gui.case.PreProcessingInformationsView import Informations, preprocessorFile
from code_saturne.model.LocalizationModel import LocalizationModel, Zone
from code_saturne.model.OutputControlModel import OutputControlModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("LocalizationView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Line edit delegate for the label
#-------------------------------------------------------------------------------

class LabelDelegate(QItemDelegate):
    """
    Use of a QLineEdit in the table.
    """
    def __init__(self, parent, mdl):
        super(LabelDelegate, self).__init__(parent)
        self.parent = parent
        self.mdl    = mdl


    def createEditor(self, parent, option, index):
        editor = QLineEdit(parent)
        rx = "[\-_A-Za-z0-9]{1," + str(LABEL_LENGTH_MAX) + "}"
        self.regExp = QRegExp(rx)
        v = RegExpValidator(editor, self.regExp)
        editor.setValidator(v)
        return editor


    def setEditorData(self, editor, index):
        editor.setAutoFillBackground(True)
        v = from_qvariant(index.model().data(index, Qt.DisplayRole), to_text_string)
        self.p_value = str(v)
        editor.setText(v)


    def setModelData(self, editor, model, index):
        if not editor.isModified():
            return

        if editor.validator().state == QValidator.Acceptable:
            p_value = str(editor.text())

            if p_value in self.mdl.getLabelsZonesList():
                default              = {}
                default['label']     = self.p_value
                default['list']      = self.mdl.getLabelsZonesList()
                default['regexp'] = self.regExp
                log.debug("setModelData-> default = %s" % default)

                from code_saturne.gui.case.VerifyExistenceLabelDialogView import VerifyExistenceLabelDialogView
                dialog = VerifyExistenceLabelDialogView(self.parent, default)
                if dialog.exec_():
                    result = dialog.get_result()
                    p_value  = result['label']
                    log.debug("setModelData-> result = %s" % result)
                else:
                    p_value = self.p_value

            model.setData(index, p_value, Qt.DisplayRole)

#-------------------------------------------------------------------------------
# QLineEdit delegate for the zone
#-------------------------------------------------------------------------------

class CodeNumberDelegate(QItemDelegate):
    def __init__(self, parent, mdl):
        super(CodeNumberDelegate, self).__init__(parent)
        self.parent = parent
        self.mdl = mdl


    def createEditor(self, parent, option, index):
        editor = QLineEdit(parent)
        validator = IntValidator(editor, min=0, max=999)
        editor.setValidator(validator)
        return editor


    def setEditorData(self, editor, index):
        editor.setAutoFillBackground(True)
        v = from_qvariant(index.model().data(index, Qt.DisplayRole), to_text_string)
        editor.setText(self.value)


    def setModelData(self, editor, model, index):
        if editor.validator().state == QValidator.Acceptable:
            value = from_qvariant(editor.text(), int)

            # Check for unicity
            if value != self.value and str(value) in self.mdl.getCodeNumbersList():
                title = self.tr("Warning")
                msg   = self.tr("Zone number can be used only once.\n"\
                            "Please give another value.")
                QMessageBox.warning(self.parent, title, msg)
                return

            model.setData(index, value, Qt.DisplayRole)


class LocalizationSelectorDelegate(QItemDelegate):
    def __init__(self, parent, mdl):
        super(LocalizationSelectorDelegate, self).__init__(parent)
        self.parent = parent
        self.mdl = mdl


    def createEditor(self, parent, option, index):

        editor = QLineEdit(parent)

        # Autocompletion for selection criteria!
        GuiLM = GuiLabelManager()
        comp_list = GuiLM.getCompleter("mesh_selection")
        comp_list += GuiLM.getCompleter("mesh_normal")

        completer = QCompleter()
        editor.setCompleter(completer)
        model = QStringListModel()
        completer.setModel(model)
        model.setStringList(comp_list)

        return editor


    def setEditorData(self, editor, index):
        # This line is used to avoid an overlay of old and new text
        editor.setAutoFillBackground(True)

        self.value = from_qvariant(index.model().data(index, Qt.DisplayRole), to_text_string)
        editor.setText(self.value)


    def setModelData(self, editor, model, index):
        value = editor.text()

        if value != self.value and str(value) in self.mdl.getLocalizationsZonesList():
            title = self.tr("Warning")
            msg   = self.tr("This localization is already used.\n"\
                            "Please give another one.")
            QMessageBox.information(self.parent, title, msg)
            return

        if str(value) != "" :
            model.setData(index, value, Qt.DisplayRole)


class DefineZonesTableModel(QStandardItemModel):
    def __init__(self, mdl, zoneType, tree=None, case=None):
        """
        """
        QStandardItemModel.__init__(self)
        self.headers = [self.tr("Label"),
                        self.tr("Zone"),
                        self.tr("Selection criteria")]
        self.setColumnCount(len(self.headers))

        self.mdl = mdl
        self.zoneType = zoneType
        self.browser = tree
        self.case = case

        self._data = []
        self._disable = []


    def data(self, index, role):
        if not index.isValid():
            return None

        if role == Qt.DisplayRole:
            row = index.row()
            col = index.column()
            return self._data[row][col]

        return None

    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled
        if (index.row(), index.column()) in self._disable:
            return Qt.ItemIsSelectable
        return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable

    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return self.headers[section]
        return None

    def setData(self, index, value, role):
        row = index.row()
        col = index.column()

        [old_label, old_code, old_local] = self._data[row]
        old_zone = self.mdl.selectZone(old_code, criterium="codeNumber")

        new_label = old_label
        new_code = old_code
        new_nature = old_zone.getNature()
        new_local = old_local

        if col == 0:
            new_label = from_qvariant(value, to_text_string)
            self._data[row][col] = new_label

        elif col == 1:
            new_code = from_qvariant(value, int)
            self._data[row][col] = new_code

        elif col == 2:
            new_local = str(from_qvariant(value, to_text_string))
            self._data[row][col] = new_local

        new_zone = Zone(self.zoneType,
                        case=self.case,
                        label=new_label,
                        codeNumber=new_code,
                        localization=new_local,
                        nature=new_nature)

        self.mdl.replaceZone(old_zone, new_zone)

        self.dataChanged.emit(index, index)
        self.browser.configureTree(self.case)
        return True

    def addItem(self, zone=None):
        """
        Add an element in the table view.
        """
        if not zone:
            zone = self.mdl.addZone(Zone(self.zoneType, case=self.case))

        line = [zone.getLabel(),
                zone.getCodeNumber(),
                zone.getLocalization()]
        self._data.append(line)
        row = self.rowCount()
        self.setRowCount(row + 1)

        # Warning: the Volume region 'all_cells' is mandatory, and can not be removed.
        if zone.getLabel() == "all_cells":
            for c in [0, 2]:
                self._disable.append((row, c))
            # self._disable.append((row, 2))
        self._disable.append((row, 1))
        self.browser.configureTree(self.case)
        return zone

    def getItem(self, row):
        return self._data[row]

    def updateItem(self):
        # update zone Id
        for id in range(0, len(self.mdl.getCodeNumbersList())):
            self._data[id][1] = id + 1

    def deleteItem(self, irow):
        del self._data[irow]
        nb_rows = self.rowCount()
        self.setRowCount(nb_rows - 1)
        self.updateItem()
        if irow < nb_rows:
            self.browser.configureTree(self.case)

    def deleteItems(self):
        for row in range(self.rowCount()):
            del self._data[0]
        self.setRowCount(0)

    def getData(self, row, column):
        return self._data[row][column]


# -------------------------------------------------------------------------------
# Main class
# -------------------------------------------------------------------------------

class LocalizationView(QWidget, Ui_LocalizationForm):
    """
    Main class
    """

    def __init__(self, zoneType, parent, case, tree=None):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_LocalizationForm.__init__(self)
        self.setupUi(self)

        self.case = case
        self.case.undoStopGlobal()

        self.zoneType = zoneType

        self.mdl = LocalizationModel(zoneType, case)
        self.case['dump_python'].append([self.mdl.__module__, zoneType, ()])

        self.browser = tree

        # Delegates
        delegateLabel = LabelDelegate(self.tableView, self.mdl)
        delegateCode = CodeNumberDelegate(self.tableView, self.mdl)
        delegateLocal = LocalizationSelectorDelegate(self.tableView, self.mdl)

        # Model for table View
        self.modelLocalization = DefineZonesTableModel(self.mdl, zoneType, tree, case)
        self.tableView.setModel(self.modelLocalization)
        self.tableView.setItemDelegateForColumn(0, delegateLabel)
        self.tableView.setItemDelegateForColumn(1, delegateCode)
        self.tableView.setItemDelegateForColumn(2, delegateLocal)
        last_section = 2

        # Populate QTableView model
        for zone in self.mdl.getZones():
            self.modelLocalization.addItem(zone)

        if QT_API == "PYQT4":
            self.tableView.verticalHeader().setResizeMode(QHeaderView.ResizeToContents)
            self.tableView.horizontalHeader().setResizeMode(QHeaderView.ResizeToContents)
            self.tableView.horizontalHeader().setResizeMode(last_section, QHeaderView.Stretch)
        elif QT_API == "PYQT5":
            self.tableView.verticalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
            self.tableView.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
            self.tableView.horizontalHeader().setSectionResizeMode(last_section, QHeaderView.Stretch)

        # Connections
        self.pushButtonNew.clicked.connect(self.slotAddZone)
        self.pushButtonDelete.clicked.connect(self.slotDeleteZone)
        self.toolButtonCreation.clicked.connect(self.slotAddFromPrePro)
        self.modelLocalization.dataChanged.connect(self.dataChanged)
        self.tableView.clicked.connect(self.slotChangeSelection)

        self.pushButtonSalome.hide()
        if case['salome']:
            self.pushButtonSalome.show()
            self.pushButtonSalome.clicked.connect(self.slotAddFromSalome)

        # Context menu
        self.tableView.setContextMenuPolicy(Qt.CustomContextMenu)
        self.tableView.customContextMenuRequested[QPoint].connect(self.slotContextMenu)

        self.case.undoStartGlobal()

    def slotChangeSelection(self):
        """
        """
        current = self.tableView.currentIndex()

    @pyqtSlot()
    def slotAddZone(self):
        """
        Insert a new item in the table view.
        """
        zone = self.modelLocalization.addItem()
        self.slotChangeSelection()


    @pyqtSlot()
    def slotDeleteZone(self):
        """
        Private Slot.
        Warning: the Volume region 'all_cells' is mandatory, therefore it can not be deleted.
        """
        lst = []
        for index in self.tableView.selectionModel().selectedRows():
            row = index.row()
            lst.append(row)

        lst.sort()
        lst.reverse()

        for row in lst:
            label = self.modelLocalization.getItem(row)[0]
            if not (label == "all_cells" and self.zoneType == 'VolumicZone'):
                # We also need to delete postprocessing meshes based on the zone
                OutputControlModel(self.case).deleteZone(label, self.zoneType)
                # Delete the zone itself
                self.mdl.deleteZone(label)
                self.modelLocalization.deleteItem(row)
        self.slotChangeSelection()
        for index in self.tableView.selectionModel().selectedRows():
            self.modelLocalization.dataChanged.emit(index, index)


    @pyqtSlot()
    def slotAddFromPrePro(self):
        """
        Research a preprocessor log to pick colors or groups of cells or faces.
        """
        if self.zoneType == 'VolumicZone':
            entity = 'cells'
        elif self.zoneType == 'BoundaryZone':
            entity = 'faces'

        file_name = preprocessorFile(self, self.case['resu_path'])

        if file_name:
            for loc in Informations(file_name, entity).getLocalizations():
                if loc not in self.mdl.getLocalizationsZonesList():
                    zone = Zone(self.zoneType, case = self.case, localization = loc)
                    self.mdl.addZone(zone)
                    self.modelLocalization.addItem(zone)


    @pyqtSlot()
    def slotContextMenu(self):
        """
        Public slot

        Create the popup menu of the Browser
        """
        fileMenu = QMenu(self.tableView)

        actionMerge = QAction(self.tr("Merge selected zones"), self.tableView)
        actionMerge.triggered.connect(self.slotMerge)
        fileMenu.addAction(actionMerge)

        fileMenu.popup(QCursor().pos())
        fileMenu.show()


    @pyqtSlot()
    def slotMerge(self):
        """
        public slot
        """
        lst = []
        for index in self.tableView.selectionModel().selectedRows():
            lst.append(index.row())

        row = lst.pop(0)
        [label, code, new_localization] = self.modelLocalization.getItem(row)
        ll = label

        for row in lst:
            [label, code, localization] = self.modelLocalization.getItem(row)
            if "all[]" not in new_localization.split(" "):
                new_localization += " or " + localization
            if localization == "all[]":
                new_localization = "all[]"

        self.modelLocalization.deleteItems()
        self.mdl.mergeZones(ll, new_localization, lst)

        # Populate QTableView model
        for zone in self.mdl.getZones():
            self.modelLocalization.addItem(zone)


    @pyqtSlot()
    def slotAddFromSalome(self):
        """
        When GUI is embeded in the Salome desktop. Add selection criteria from
        graphical selection in the VTK viwver, or in the ObjectBrowser.
        """
        if self.case['salome']:
            from code_saturne.gui.case.SalomeHandler import BoundaryGroup, VolumeGroup

            log.debug("slotAddFromSalome: zoneType -> %s" % self.zoneType)
            if self.zoneType == 'VolumicZone':
                loc = VolumeGroup()
            elif self.zoneType == 'BoundaryZone':
                loc = BoundaryGroup()

            log.debug("slotAddFromSalome: selection criteria -> %s" % loc)
            if loc not in self.mdl.getLocalizationsZonesList() and loc != "" :
                zone = Zone(self.zoneType, localization = loc, case = self.case)
                self.mdl.addZone(zone)
                self.modelLocalization.addItem(zone)


    def dataChanged(self, topLeft, bottomRight):
        for row in range(topLeft.row(), bottomRight.row()+1):
            self.tableView.resizeRowToContents(row)
        for col in range(topLeft.column(), bottomRight.column()+1):
            self.tableView.resizeColumnToContents(col)


#-------------------------------------------------------------------------------
# Define Volumic regions class
#-------------------------------------------------------------------------------

class VolumeLocalizationView(LocalizationView):
    """
    Define volume regions class.
    """
    def __init__(self, parent=None, case=None, tree = None, hide_all=False):
        """
        Constructor.
        """
        self.case = case
        self.browser = tree

        LocalizationView.__init__(self, "VolumicZone", parent, self.case, tree)

        title = self.tr("Definition of volume regions")
        self.groupBoxLocalization.setTitle(title)

        if hide_all:
            self.groupBoxLocalization.hide()

        self.browser.configureTree(self.case)

#-------------------------------------------------------------------------------
# Define Boundary regions class
#-------------------------------------------------------------------------------

class BoundaryLocalizationView(LocalizationView):
    """
    Define boundary regions class.
    """
    def __init__(self, parent=None, case=None, tree = None):
        """
        Constructor.
        """
        self.case = case
        self.browser = tree

        LocalizationView.__init__(self, "BoundaryZone", parent, self.case, tree)

        title = self.tr("Boundary regions definition")
        self.groupBoxLocalization.setTitle(title)


#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
