# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2023 EDF S.A.
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
This module contains the following classes and function:
- StandardItemModelMeshes
- SolutionDomainView
"""

#-------------------------------------------------------------------------------
# Standard modules
#-------------------------------------------------------------------------------

import os, sys, logging
import configparser

#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------

from code_saturne.gui.base.QtCore    import *
from code_saturne.gui.base.QtGui     import *
from code_saturne.gui.base.QtWidgets import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.model.Common import GuiParam
from code_saturne.gui.base.QtPage import ComboModel, DoubleValidator, RegExpValidator, IntValidator, BasicTableModel, LabelDelegate
from code_saturne.gui.base.QtPage import from_qvariant, to_text_string
from code_saturne.gui.case.TimeTablesForm import Ui_TimeTablesForm
from code_saturne.model.SolutionDomainModel import RelOrAbsPath, MeshModel, SolutionDomainModel
from code_saturne.model.TimeTablesModel import TimeTablesModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("TimeTablesView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Label delegate for 'Name' in Meshes table
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

class StandardItemModelTables(QStandardItemModel):

    def __init__(self, mdl):
        """
        Constructor
        """
        QStandardItemModel.__init__(self)
        self.mdl = mdl
        self.dataTables = []

        self.populateModel()

        self.headers = [self.tr("Table name"),
                        self.tr("Separator"),
                        self.tr("File name")]

        self.tooltip = [self.tr("Name of the table"),
                        self.tr("Character separating data within a line"),
                        self.tr("Table file name")]

        self.setColumnCount(len(self.headers))



    def populateModel(self):

        for i in range(self.mdl.getNumberOfTables()):
            _d = [self.mdl.getTableName(i),
                  self.mdl.getTableDelimiter(i),
                  self.mdl.getTableFileName(i)]

            self.dataTables.append(_d)
            log.debug("populateModel -> dataTables = %s" % _d)
            row = self.rowCount()
            self.setRowCount(row + 1)


    def data(self, index, role):
        """
        data method
        """
        if not index.isValid():
            return None

        col = index.column()

        if role == Qt.ToolTipRole:
            return self.tooltip[col]

        elif role == Qt.DisplayRole:
            d = self.dataTables[index.row()][col]
            if d:
                return d
            else:
                return None

        return None


    def flags(self, index):
        """
        Return Qt flags
        """
        if not index.isValid():
            return Qt.ItemIsEnabled

        if index.column() == 2:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable
        else:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable


    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return self.headers[section]
        else:
            return None

    def setData(self, index, value, role=None):
        row = index.row()
        col = index.column()

        v = from_qvariant(value, to_text_string)
        if col == 0:
            self.dataTables[row][col] = v
            if v is not None:
                self.mdl.setTableName(row, v)
        elif col == 1:
            self.dataTables[row][col] = v
            if v is not None:
                self.mdl.setTableDelimiter(row, v)

        self.dataChanged.emit(index, index)
        return True

    def addRow(self, table_file_name):
        """
        """
        _default_name = os.path.basename(table_file_name).split('.')[0]
        self.dataTables.append([_default_name,
                                self.mdl.defaultValues('delimiter'),
                                table_file_name])

        row = self.rowCount()
        self.setRowCount(row+1)

        self.mdl.addTable()
        self.mdl.setTableName(row, _default_name)
        self.mdl.setTableDelimiter(row, self.mdl.defaultValues('delimiter'))
        self.mdl.setTableFileName(row, table_file_name)


    def deleteRow(self, row):
        """
        Delete a row in the model
        """
        self.setRowCount(self.rowCount() - 1)
        del self.dataTables[row]


#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class TimeTablesView(QWidget, Ui_TimeTablesForm):
    """
    """
    def __init__(self, parent, case, stbar, tree):
        """
        Constructor
        """
        QWidget.__init__(self, parent)
        Ui_TimeTablesForm.__init__(self)
        self.setupUi(self)

        self.root = parent
        self.case = case
        self.stbar = stbar
        self.tree = tree

        self.mdl = TimeTablesModel(self.case)

        self.case.undoStopGlobal()

        self.table_id = -1

        # Set parameters for
        self.modelColumnsImport = ComboModel(self.comboBoxCols, 2, 1)
        for opt in ('all', 'subset'):
            self.modelColumnsImport.addItem(self.tr(opt, opt))

        self.modelHeadersImport = ComboModel(self.comboBoxHeaders, 2, 1)
        for opt in ('import', 'user'):
            self.modelHeadersImport.addItem(self.tr(opt, opt))


        self.groupBoxTableParameters.hide()

        self.tableModelTimeTables = StandardItemModelTables(self.mdl)


        self.tableViewTimeTables.setModel(self.tableModelTimeTables)
        self.tableViewTimeTables.resizeRowsToContents()
        self.tableViewTimeTables.setAlternatingRowColors(True)
        self.tableViewTimeTables.setSelectionBehavior(QAbstractItemView.SelectItems)
        self.tableViewTimeTables.setSelectionMode(QAbstractItemView.ExtendedSelection)

        nameDelegate = LabelDelegate(self.tableViewTimeTables,
                                     xml_model=self.mdl)
        self.tableViewTimeTables.setItemDelegateForColumn(0, nameDelegate)

        delimDelegate = LabelDelegate(self.tableViewTimeTables,
                                      xml_model=self.mdl)
        self.tableViewTimeTables.setItemDelegateForColumn(1, delimDelegate)

        fpathDelegate = LabelDelegate(self.tableViewTimeTables,
                                      xml_model=self.mdl)
        self.tableViewTimeTables.setItemDelegateForColumn(2, fpathDelegate)

        self.tableViewTimeTables.resizeColumnToContents(2)

        if QT_API == "PYQT4":
            self.tableViewTimeTables.verticalHeader().setResizeMode(QHeaderView.ResizeToContents)
            self.tableViewTimeTables.horizontalHeader().setResizeMode(QHeaderView.ResizeToContents)
            self.tableViewTimeTables.horizontalHeader().setResizeMode(2, QHeaderView.Stretch)
        elif QT_API == "PYQT5":
            self.tableViewTimeTables.verticalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
            self.tableViewTimeTables.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
            self.tableViewTimeTables.horizontalHeader().setSectionResizeMode(2, QHeaderView.Stretch)

        #Connections
        self.pushButtonAddTimeTable.clicked.connect(self.slotSearchTable)
        self.pushButtonDeleteTimeTable.clicked.connect(self.slotDeleteTable)
        self.tableViewTimeTables.clicked.connect(self.slotSelectTable)

        self.spinBoxRows.valueChanged[int].connect(self.slotRowsToSkip)
        self.comboBoxCols.activated[str].connect(self.slotColImportMode)
        self.lineEditColumnsList.textChanged[str].connect(self.slotColToImport)

        self.comboBoxHeaders.activated[str].connect(self.slotHeadImportMode)
        self.lineEditHeadersList.textChanged[str].connect(self.slotHeadToImport)

        self.pushButtonImportheaders.clicked.connect(self.slotImportHeadersFromFile)
        self.spinBoxHeadersRow.valueChanged[int].connect(self.slotHeadLine)

    def _update_page_view(self):
        """
        Update state of widgets
        """
        if self.table_id < 0:
            return

        self.groupBoxTableParameters.show()

        # Columns
        import_all_cols = \
                self.mdl.getTableProperty(self.table_id, 'cols2import') == 'all'

        self.lineEditColumnsList.setEnabled(not import_all_cols)

        # Headers
        import_headers = \
                self.mdl.getTableProperty(self.table_id, 'headers_def') == 'import'

        self.lineEditHeadersList.setEnabled(not import_headers)
        self.spinBoxHeadersRow.setEnabled(import_headers)

        self.pushButtonImportheaders.setEnabled(import_headers)


    @pyqtSlot("QModelIndex")
    def slotSelectTable(self, index):
        """
        Action when user clicks on the table
        """
        cindex = self.tableViewTimeTables.currentIndex()
        if cindex != (-1, -1):
            self.table_id = cindex.row()
            self.spinBoxRows.setValue(int(self._get_mdl_property('skip_rows')))
            self.comboBoxCols.setCurrentText(self._get_mdl_property('cols2import'))
            self.lineEditColumnsList.setText(self._get_mdl_property('col_ids'))
            self.comboBoxHeaders.setCurrentText(self._get_mdl_property('headers_def'))
            self.lineEditHeadersList.setText(self._get_mdl_property('headers_list'))
            self.spinBoxHeadersRow.setValue(int(self._get_mdl_property('headers_line')))

            self._update_page_view()


    @pyqtSlot()
    def slotSearchTable(self):
        """
        Slot for file selection
        """
        msg = self.tr("Select a time table file.")
        self.stbar.showMessage(msg, 2000)

        table_files = self.selectTableFile()

        for _f in table_files:
            self.tableModelTimeTables.addRow(_f)


    @pyqtSlot()
    def slotDeleteTable(self):
        """
        Delete an entry in the table
        """
        current = self.tableViewTimeTables.currentIndex()
        idx = current.row()
        self.tableModelTimeTables.deleteRow(idx)
        self.mdl.deleteTable(idx)


    @pyqtSlot(int)
    def slotRowsToSkip(self, val):
        """
        Set number of rows to skip when importing table
        """
        self.mdl.setTableProperty(self.table_id, 'skip_rows', str(val))

    @pyqtSlot(str)
    def slotColImportMode(self, text):
        """
        Set import mode of data columns
        """
        self.mdl.setTableProperty(self.table_id, 'cols2import', str(text))
        self._update_page_view()

    @pyqtSlot(str)
    def slotColToImport(self, text):
        """
        Set selection criteria for columns (if user defined import mode)
        """
        self.mdl.setTableProperty(self.table_id, 'col_ids', str(text))

    @pyqtSlot(str)
    def slotHeadImportMode(self, text):
        """
        Set import mode for headers (data labels)
        """
        self.mdl.setTableProperty(self.table_id, 'headers_def', str(text))

        self.lineEditHeadersList.setReadOnly(text == 'file')

        self._update_page_view()

    @pyqtSlot(str)
    def slotHeadToImport(self, text):
        """
        Set headers value
        """
        self.mdl.setTableProperty(self.table_id, 'headers_list', str(text))


    @pyqtSlot()
    def slotImportHeadersFromFile(self):
        """
        Import headers fro a data file
        """
        fname = self.mdl.getTableFileName(self.table_id)
        d = self.mdl.getTableDelimiter(self.table_id)

        hl = self._get_headers_from_file(fname, d)

        if hl != "":
            self.mdl.setTableProperty(self.table_id, 'headers_list', hl)

        self.lineEditHeadersList.setText(hl)
        return


    @pyqtSlot(int)
    def slotHeadLine(self, val):
        """
        Define line from which to import headers
        """
        self.mdl.setTableProperty(self.table_id, 'headers_line', int(val))


    def selectTableFile(self):
        """
        Open a File dialog in order to select a time table.
        """
        _files = []

        title = self.tr("Select input file for time table.")

        default = os.path.split(self.case['case_path'])[0]

        if hasattr(QFileDialog, 'ReadOnly'):
            options = QFileDialog.DontUseNativeDialog | QFileDialog.ReadOnly
        else:
            options = QFileDialog.DontUseNativeDialog

        filetypes = self.tr('csv file (*.csv);;text file (*.txt);;data file (*.dat);;All files (*)')

        dialog = QFileDialog()
        dialog.setWindowTitle(title)
        dialog.setDirectory(default)
        dialog.setNameFilter(filetypes)

        if hasattr(dialog, 'setOptions'):
            dialog.setOptions(options)
        dialog.setFileMode(QFileDialog.ExistingFiles)

        if dialog.exec_() == 1:
            s = dialog.selectedFiles()
            count = len(s)
            for i in range(count):
                el = str(s[0])
                s = s[1:]
                _files.append(el)

        return _files


    def _get_headers_from_file(self, file_name, delim):
        """
        Get headers from a file
        """
        id_head = int(self.mdl.getTableProperty(self.table_id, 'headers_line')) - 1

        retval = ""
        with open(file_name) as fp:
            lines = fp.readlines()
            if id_head > len(lines):
                return retval

            if len(delim.strip()) > 0:
                headers = [elt.strip() for elt in lines[id_head].split(delim)]
            else:
                headers = [elt.strip() for elt in lines[id_head].split()]

            for i, h in enumerate(headers):
                if i != 0:
                    retval += ","
                retval += str(h)

        return retval

    def _get_mdl_property(self, ppty):
        """
        """
        return self.mdl.getTableProperty(self.table_id, ppty)
