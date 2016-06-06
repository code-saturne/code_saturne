# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2016 EDF S.A.
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
- LabelDelegate
- CodeNumberDelegate
- BoundaryNatureDelegate
- StandardItemVolumeNature
- FlagBox
- VolumeNatureDelegate
- LocalizationSelectorDelegate
- StandardItemModelLocalization
- LocalizationView
- VolumeLocalizationView
- BoundaryLocalizationView
"""

#-------------------------------------------------------------------------------
# Standard modules
#-------------------------------------------------------------------------------

import string, logging

#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------

from code_saturne.Base.QtCore    import *
from code_saturne.Base.QtGui     import *
from code_saturne.Base.QtWidgets import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.Base.Toolbox import GuiParam
from code_saturne.Base.Common import LABEL_LENGTH_MAX
from code_saturne.Base.QtPage import IntValidator, RegExpValidator
from code_saturne.Base.QtPage import to_qvariant, from_qvariant, to_text_string
from code_saturne.Pages.LocalizationForm import Ui_LocalizationForm
from code_saturne.Pages.PreProcessingInformationsView import Informations, preprocessorFile
from code_saturne.Pages.LocalizationModel import LocalizationModel, Zone

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

                from code_saturne.Pages.VerifyExistenceLabelDialogView import VerifyExistenceLabelDialogView
                dialog = VerifyExistenceLabelDialogView(self.parent, default)
                if dialog.exec_():
                    result = dialog.get_result()
                    p_value  = result['label']
                    log.debug("setModelData-> result = %s" % result)
                else:
                    p_value = self.p_value

            model.setData(index, to_qvariant(p_value), Qt.DisplayRole)

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

            model.setData(index, to_qvariant(value), Qt.DisplayRole)

#-------------------------------------------------------------------------------
# QComboBox delegate for the boundary nature
#-------------------------------------------------------------------------------

class BoundaryNatureDelegate(QItemDelegate):
    """
    Use of a combo box in the table.
    """
    def __init__(self, parent, dicoM2V):
        super(BoundaryNatureDelegate, self).__init__(parent)
        self.parent   = parent
        self.dicoM2V = dicoM2V

        self.dicoV2M = {}
        for k, v in list(self.dicoM2V.items()):
            self.dicoV2M[v] = k


    def createEditor(self, parent, option, index):
        editor = QComboBox(parent)
        for k in list(self.dicoV2M.keys()):
            editor.addItem(k)
        editor.installEventFilter(self)
        editor.setMinimumWidth(120)
        return editor


    def setEditorData(self, comboBox, index):
        row = index.row()
        col = index.column()
        str_model = index.model().getData(row, col)
        idx = list(self.dicoM2V.keys()).index(str_model)
        comboBox.setCurrentIndex(idx)


    def setModelData(self, comboBox, model, index):
        txt   = str(comboBox.currentText())
        value = self.dicoV2M[txt]
        selectionModel = self.parent.selectionModel()
        for idx in selectionModel.selectedIndexes():
            if idx.column() == index.column():
                model.setData(idx, to_qvariant(value), Qt.DisplayRole)


    def tr(self, text):
        return text

#-------------------------------------------------------------------------------
# StandarItemModel class for QComboBox in nature delegate (Volumic)
#-------------------------------------------------------------------------------

class StandardItemVolumeNature(QStandardItemModel):
    """
    QStandardItemModel associated with the QComboBox
    in the delegate (validator) for the nature in the case of volumic zone.

    Attr. natureList is a dictionnary between model string and view string.
    Attr. dicoNature holds value 'on' or 'off' if the item is checked or not.
    """

    def __init__(self, dicoM2V, case, data):
        """
        """
        QStandardItemModel.__init__(self)

        # Dictionary item to label for display
        self.dicoM2V = dicoM2V
        # Dictionary item to status (on/off)
        self.dicoNature = data

        #self.keys = self.natureList.keys() # ordered tuple
        self.keys = Zone('VolumicZone', case = case).getNatureList()

        self.setColumnCount(1)

        row = 0
        for key in self.keys:
            self.setItem(row, QStandardItem(str(self.dicoM2V[key])))
            row += 1
        self.setRowCount(row)


    def data(self, index, role):
        if not index.isValid():
            return to_qvariant()

        if role == Qt.EditRole:
            key = self.keys[index.row()]
            return to_qvariant(self.dicoM2V[key])

        elif role == Qt.DisplayRole:
            key = self.keys[index.row()]
            return to_qvariant(self.dicoM2V[key])

        elif role == Qt.CheckStateRole:
            key = self.keys[index.row()]
            value = self.dicoNature[key]
            if value == 'on':
                return to_qvariant(Qt.Checked)
            else:
                return to_qvariant(Qt.Unchecked)

        return to_qvariant()


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled

        return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsUserCheckable


    def setData(self, index, value, role=None):
        if not index.isValid():
            return

        if role == Qt.EditRole:
            pass

        elif role == Qt.DisplayRole:
            key = self.keys[index.row()]
            self.dicoNature[key] = str(value)

        elif role == Qt.CheckStateRole:
            state = from_qvariant(value, int)
            key = self.keys[index.row()]
            if state == Qt.Unchecked:
                self.dicoNature[key] = "off"
            else:
                self.dicoNature[key] = "on"

        id1 = self.index(0, 0)
        id2 = self.index(self.rowCount(), 0)
        self.dataChanged.emit(id1, id2)
        return True


    def getChecked(self):
        s = []
        for k, v in list(self.dicoNature.items()):
            if v == "on":
                s.append(k)
        return ";".join(s)

#-------------------------------------------------------------------------------
# FlagBox: new QComboBox to construct Delegate for the volume nature
#-------------------------------------------------------------------------------

class FlagBox(QComboBox):
    def __init__(self, parent):
        QComboBox.__init__(self, parent)

        opt = QStyleOptionComboBox()
        opt.initFrom(self)
        opt.editable = self.isEditable()
        if (self.style().styleHint(QStyle.SH_ComboBox_Popup, opt)):
            self.setItemDelegate(QItemDelegate(self))

        self.activated[int].connect(self.slotActivated)


    @pyqtSlot(int)
    def slotActivated(self, index):
        value = self.itemData(index, Qt.CheckStateRole)
        state = from_qvariant(value, int)
        if state == Qt.Unchecked:
            s = Qt.Checked
        else:
            s = Qt.Unchecked
        self.setItemData(index, to_qvariant(s), Qt.CheckStateRole)


#-------------------------------------------------------------------------------
# Delegate for the volume nature
#-------------------------------------------------------------------------------

class VolumeNatureDelegate(QItemDelegate):
    """
    Use of a combo box in the table.
    """
    def __init__(self, parent, case):
        super(VolumeNatureDelegate, self).__init__(parent)
        self.case = case


    def createEditor(self, parent, option, index):
        editor = FlagBox(parent)
        editor.setEditable(False)
        editor.setMinimumWidth(160)
        editor.installEventFilter(self)
        return editor


    def setEditorData(self, editor, index):
        data = index.model().getData(index.row(), index.column())
        dicoM2V = index.model().dicoM2V

        self.flagbox_model = StandardItemVolumeNature(dicoM2V, self.case, data)
        editor.setModel(self.flagbox_model)


    def setModelData(self, editor, model, index):
        value = self.flagbox_model.getChecked()
        model.setData(index, to_qvariant(value), Qt.DisplayRole)

#-------------------------------------------------------------------------------
# QLineEdit delegate for localization
#-------------------------------------------------------------------------------

class LocalizationSelectorDelegate(QItemDelegate):
    def __init__(self, parent, mdl):
        super(LocalizationSelectorDelegate, self).__init__(parent)
        self.parent = parent
        self.mdl = mdl


    def createEditor(self, parent, option, index):
        editor = QLineEdit(parent)
        return editor


    def setEditorData(self, editor, index):
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
            model.setData(index, to_qvariant(value), Qt.DisplayRole)

#-------------------------------------------------------------------------------
# StandarItemModel class
#-------------------------------------------------------------------------------

class StandardItemModelLocalization(QStandardItemModel):
    def __init__(self, mdl, zoneType, dicoM2V, tree = None, case = None):
        """
        """
        QStandardItemModel.__init__(self)
        self.headers = [self.tr("Label"),
                        self.tr("Zone"),
                        self.tr("Nature"),
                        self.tr("Selection criteria")]
        self.setColumnCount(len(self.headers))

        self.mdl      = mdl
        self.zoneType = zoneType
        self.dicoM2V  = dicoM2V
        self.browser = tree
        self.case = case

        self._data = []
        self._disable = []


    def data(self, index, role):
        if not index.isValid():
            return to_qvariant()

        if role == Qt.DisplayRole:
            row = index.row()
            col = index.column()

            if col in [0, 1, 3]:
                return to_qvariant(self._data[row][col])

            elif col == 2:
                if self.zoneType == "VolumicZone":
                    data = self._data[row][col]
                    item = "\n".join([self.dicoM2V[key] for key in list(self.dicoM2V.keys()) if data[key] == "on"])

                    return to_qvariant(item)

                elif self.zoneType == "BoundaryZone":
                    key = self._data[row][col]
                    return to_qvariant(self.dicoM2V[key])

        return to_qvariant()


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled
        if (index.row(), index.column()) in self._disable:
            return Qt.ItemIsSelectable
        return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable


    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return to_qvariant(self.headers[section])
        return to_qvariant()


    def setData(self, index, value, role):
        row = index.row()
        col = index.column()

        [old_label, old_code, old_nature, old_local] = self._data[row]

        old_zone = Zone(self.zoneType,
                        case         = self.case,
                        label        = old_label,
                        codeNumber   = old_code,
                        localization = old_local,
                        nature       = old_nature)

        new_label  = old_label
        new_code   = old_code
        new_nature = old_nature
        new_local  = old_local

        if col == 0:
            new_label = from_qvariant(value, to_text_string)
            self._data[row][col] = new_label

        elif col == 1:
            new_code = from_qvariant(value, int)
            self._data[row][col] = new_code

        elif col == 2:

            if self.zoneType == "VolumicZone":
                # We modify the dictionary here
                nature_list = str(from_qvariant(value, to_text_string)).split(";")

                for key in list(self._data[row][col].keys()):
                    if key in nature_list:
                        self._data[row][col][key] = "on"
                    else:
                        self._data[row][col][key] = "off"

                    new_nature = self._data[row][col].copy()

            elif self.zoneType == "BoundaryZone":
                new_nature = str(from_qvariant(value, to_text_string))
                self._data[row][col] = new_nature

        elif col == 3:
            new_local = str(from_qvariant(value, to_text_string))
            self._data[row][col] = new_local

        new_zone = Zone(self.zoneType,
                        case         = self.case,
                        label        = new_label,
                        codeNumber   = new_code,
                        localization = new_local,
                        nature       = new_nature)

        self.mdl.replaceZone(old_zone, new_zone)

        self.dataChanged.emit(index, index)
        self.browser.configureTree(self.case)
        return True


    def addItem(self, zone=None):
        """
        Add an element in the table view.
        """
        if not zone:
            zone = self.mdl.addZone(Zone(self.zoneType, case = self.case))

        line = [zone.getLabel(),
                zone.getCodeNumber(),
                zone.getNature(),
                zone.getLocalization()]
        self._data.append(line)
        row = self.rowCount()
        self.setRowCount(row+1)

        # Warning: the Volume region 'all_cells' is mandatory, and can not be modified.
        if self.zoneType == "VolumicZone":
            if zone.getLabel() == "all_cells":
                for c in range(self.columnCount()):
                    self._disable.append((row, c))
        self._disable.append((row, 1))
        self.browser.configureTree(self.case)


    def getItem(self, row):
        return self._data[row]


    def updateItem(self):
        # update zone Id
        for id in range(0, len(self.mdl.getCodeNumbersList())):
            self._data[id][1] = id + 1


    def deleteItem(self, irow):
        del self._data[irow]
        nb_rows = self.rowCount()
        self.setRowCount(nb_rows-1)
        self.updateItem()


    def deleteItems(self):
        for row in range(self.rowCount()):
            del self._data[0]
        self.setRowCount(0)


    def getData(self, row, column):
        return self._data[row][column]

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class LocalizationView(QWidget, Ui_LocalizationForm):
    """
    Main class
    """
    def __init__(self, zoneType, parent, case, dicoM2V, tree = None):
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

        # Model for table View
        self.modelLocalization = StandardItemModelLocalization(self.mdl, zoneType, dicoM2V, tree, case)
        self.tableView.setModel(self.modelLocalization)

        # Delegates
        delegateLabel = LabelDelegate(self.tableView, self.mdl)
        delegateCode  = CodeNumberDelegate(self.tableView, self.mdl)
        delegateLocal = LocalizationSelectorDelegate(self.tableView, self.mdl)

        self.tableView.setItemDelegateForColumn(0, delegateLabel)
        self.tableView.setItemDelegateForColumn(1, delegateCode)
        self.tableView.setItemDelegateForColumn(3, delegateLocal)

        # Populate QTableView model
        for zone in self.mdl.getZones():
            self.modelLocalization.addItem(zone)

        if QT_API == "PYQT4":
            self.tableView.verticalHeader().setResizeMode(QHeaderView.ResizeToContents)
            self.tableView.horizontalHeader().setResizeMode(QHeaderView.ResizeToContents)
            self.tableView.horizontalHeader().setResizeMode(3, QHeaderView.Stretch)
        elif QT_API == "PYQT5":
            self.tableView.verticalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
            self.tableView.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
            self.tableView.horizontalHeader().setSectionResizeMode(3, QHeaderView.Stretch)

        # Connections
        self.pushButtonNew.clicked.connect(self.slotAddZone)
        self.pushButtonDelete.clicked.connect(self.slotDeleteZone)
        self.toolButtonCreation.clicked.connect(self.slotAddFromPrePro)
        self.modelLocalization.dataChanged.connect(self.dataChanged)

        self.pushButtonSalome.hide()
        if case['salome']:
            self.pushButtonSalome.show()
            self.pushButtonSalome.clicked.connect(self.slotAddFromSalome)

        # Context menu
        self.tableView.setContextMenuPolicy(Qt.CustomContextMenu)
        self.tableView.customContextMenuRequested[QPoint].connect(self.slotContextMenu)

        self.case.undoStartGlobal()


    @pyqtSlot()
    def slotAddZone(self):
        """
        Insert a new item in the table view.
        """
        self.modelLocalization.addItem()


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
            [label, codeNumber, nature, localization] = self.modelLocalization.getItem(row)
            if not (label == "all_cells" and self.zoneType == 'VolumicZone'):
                self.mdl.deleteZone(label)
                self.modelLocalization.deleteItem(row)


    @pyqtSlot()
    def slotAddFromPrePro(self):
        """
        Research a Preprocessor listing to pick colors or groups of cells or faces.
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

        if self.zoneType == 'BoundaryZone':
            fileMenu.addSeparator()

            self.actionInlet = QAction(self.tr("Select all inlets"), self.tableView)
            self.actionInlet.triggered.connect(self.slotSelectBoudaries)
            fileMenu.addAction(self.actionInlet)

            self.actionOutlet = QAction(self.tr("Select all outlets"), self.tableView)
            self.actionOutlet.triggered.connect(self.slotSelectBoudaries)
            fileMenu.addAction(self.actionOutlet)

            self.actionWall = QAction(self.tr("Select all walls"), self.tableView)
            self.actionWall.triggered.connect(self.slotSelectBoudaries)
            fileMenu.addAction(self.actionWall)

            self.actionSymmetry = QAction(self.tr("Select all symmetries"), self.tableView)
            self.actionSymmetry.triggered.connect(self.slotSelectBoudaries)
            fileMenu.addAction(self.actionSymmetry)

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
        [label, code, nature, new_localization] = self.modelLocalization.getItem(row)
        ll = label

        for row in lst:
            [label, code, nature, localization] = self.modelLocalization.getItem(row)
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
            from code_saturne.Pages.SalomeHandler import BoundaryGroup, VolumeGroup

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


    def tr(self, text):
        """
        Translation.
        """
        return text

#-------------------------------------------------------------------------------
# Define Volumic regions class
#-------------------------------------------------------------------------------

class VolumeLocalizationView(LocalizationView):
    """
    Define volume regions class.
    """
    def __init__(self, parent=None, case=None, tree = None):
        """
        Constructor.
        """
        self.case = case
        self.browser = tree
        dicoM2V = Zone("VolumicZone", self.case).getModel2ViewDictionary()

        LocalizationView.__init__(self, "VolumicZone", parent, self.case, dicoM2V, tree)

        title = self.tr("Definition of volume regions")
        self.groupBoxLocalization.setTitle(title)

        # Delegates
        delegateNature = VolumeNatureDelegate(self.tableView, self.case)
        self.tableView.setItemDelegateForColumn(2, delegateNature)
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
        dicoM2V = Zone("BoundaryZone", case = self.case).getModel2ViewDictionary()
        self.browser = tree

        LocalizationView.__init__(self, "BoundaryZone", parent, self.case, dicoM2V, tree)

        title = self.tr("Definition of boundary regions")
        self.groupBoxLocalization.setTitle(title)

        # Delegates
        delegateNature = BoundaryNatureDelegate(self.tableView, dicoM2V)
        self.tableView.setItemDelegateForColumn(2, delegateNature)


    @pyqtSlot()
    def slotSelectBoudaries(self):
        """
        Public slot.

        Warning: works only if the selection mode of the view is set to MultiSelection.
        """
        previous_selecion_mode = self.tableView.selectionMode()
        self.tableView.setSelectionMode(QAbstractItemView.MultiSelection)
        self.tableView.clearSelection()

        if self.sender() == self.actionInlet:
            select  = "inlet"
        elif self.sender() == self.actionOutlet:
            select  = "outlet"
        elif self.sender() == self.actionWall:
            select  = "wall"
        elif self.sender() == self.actionSymmetry:
            select  = "symmetry"

        for row in range(self.modelLocalization.rowCount()):
            [label, code, nature, localization] = self.modelLocalization.getItem(row)
            if nature == select:
                self.tableView.selectRow(row)

        self.tableView.setSelectionMode(previous_selecion_mode)

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
