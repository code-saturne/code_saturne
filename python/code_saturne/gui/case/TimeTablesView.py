# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2026 EDF S.A.
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
from code_saturne.gui.case.TimeTablesForm import Ui_TimeTablesForm
from code_saturne.gui.case.TimeTablePlotDialog import TimeTablePlotDialog
from code_saturne.model.SolutionDomainModel import RelOrAbsPath, MeshModel, SolutionDomainModel
from code_saturne.model.TimeTablesModel import TimeTablesModel

if QT_API == "PYQT6":
    from code_saturne.gui.case import resources_pages_rc

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

        if role == Qt.ItemDataRole.ToolTipRole:
            return self.tooltip[col]

        elif role == Qt.ItemDataRole.DisplayRole:
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
            return Qt.ItemFlag.ItemIsEnabled

        if index.column() == 2:
            return Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsSelectable
        else:
            return Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsSelectable | Qt.ItemFlag.ItemIsEditable


    def headerData(self, section, orientation, role):
        if orientation == Qt.Orientation.Horizontal and role == Qt.ItemDataRole.DisplayRole:
            return self.headers[section]
        else:
            return None

    def setData(self, index, value, role=None):
        row = index.row()
        col = index.column()

        v = str(value)
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

        self.modelTimeOffset = ComboModel(self.comboBoxTimeOffset, 2, 1)
        for opt in ('no', 'yes'):
            self.modelTimeOffset.addItem(self.tr(opt, opt))


        self.groupBoxTableParameters.hide()

        self.tableModelTimeTables = StandardItemModelTables(self.mdl)


        self.tableViewTimeTables.setModel(self.tableModelTimeTables)
        self.tableViewTimeTables.resizeRowsToContents()
        self.tableViewTimeTables.setAlternatingRowColors(True)
        self.tableViewTimeTables.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectItems)
        self.tableViewTimeTables.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)

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

        self.tableViewTimeTables.verticalHeader().setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
        self.tableViewTimeTables.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
        self.tableViewTimeTables.horizontalHeader().setSectionResizeMode(2, QHeaderView.ResizeMode.Stretch)

        #Connections
        self.pushButtonAddTimeTable.clicked.connect(self.slotSearchTable)
        self.pushButtonDeleteTimeTable.clicked.connect(self.slotDeleteTable)
        self.pushButtonPlotTimeTable.clicked.connect(self.slotPlotTable)
        self.pushButtonPlotTimeTable.setEnabled(False)
        self.tableViewTimeTables.clicked.connect(self.slotSelectTable)

        self.spinBoxRows.valueChanged[int].connect(self.slotRowsToSkip)
        self.comboBoxCols.activated[int].connect(self.slotColImportMode)
        self.lineEditColumnsList.textChanged[str].connect(self.slotColToImport)

        self.comboBoxHeaders.activated[int].connect(self.slotHeadImportMode)
        self.lineEditHeadersList.textChanged[str].connect(self.slotHeadToImport)

        self.pushButtonImportheaders.clicked.connect(self.slotImportHeadersFromFile)
        self.spinBoxHeadersRow.valueChanged[int].connect(self.slotHeadLine)

        self.comboBoxTimeOffset.activated[int].connect(self.slotTimeOffsetMode)
        self.lineEditTimeOffset.textChanged[str].connect(self.slotTimeOffsetVal)
        self.lineEditTimeOffset.setValidator(DoubleValidator(self.lineEditTimeOffset))

    def _update_page_view(self):
        """
        Update state of widgets
        """
        if self.table_id < 0:
            self.pushButtonPlotTimeTable.setEnabled(False)
            return

        self.pushButtonPlotTimeTable.setEnabled(True)
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

        # Time offset
        use_offset = \
                self.mdl.getTableProperty(self.table_id, 'time_offset_choice') != 'no'

        self.lineEditTimeOffset.setEnabled(use_offset)


    @Slot("QModelIndex")
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
            self.comboBoxTimeOffset.setCurrentText(self._get_mdl_property('time_offset_choice'))
            self.lineEditTimeOffset.setText(str(self._get_mdl_property('time_offset_value')))

            self._update_page_view()


    @Slot()
    def slotSearchTable(self):
        """
        Slot for file selection
        """
        msg = self.tr("Select a time table file.")
        self.stbar.showMessage(msg, 2000)

        table_files = self.selectTableFile()

        _datap = os.path.abspath(os.path.join(self.case['case_path'], 'DATA'))
        for _f in table_files:
            if os.path.abspath(os.path.dirname(_f)) == _datap:
                self.tableModelTimeTables.addRow(os.path.basename(_f))
            else:
                self.tableModelTimeTables.addRow(_f)


    @Slot()
    def slotDeleteTable(self):
        """
        Delete an entry in the table
        """
        current = self.tableViewTimeTables.currentIndex()
        idx = current.row()
        if idx >= 0:
            self.tableModelTimeTables.deleteRow(idx)
            self.mdl.deleteTable(idx)
            self.table_id = -1
            self.groupBoxTableParameters.hide()
            self.pushButtonPlotTimeTable.setEnabled(False)


    @Slot(int)
    def slotRowsToSkip(self, val):
        """
        Set number of rows to skip when importing table
        """
        self.mdl.setTableProperty(self.table_id, 'skip_rows', str(val))

    @Slot(int)
    def slotColImportMode(self, idx):
        """
        Set import mode of data columns
        """
        text = self.comboBoxCols.currentText()
        self.mdl.setTableProperty(self.table_id, 'cols2import', str(text))
        self._update_page_view()

    @Slot(str)
    def slotColToImport(self, text):
        """
        Set selection criteria for columns (if user defined import mode)
        """
        self.mdl.setTableProperty(self.table_id, 'col_ids', str(text))

    @Slot(int)
    def slotHeadImportMode(self, idx):
        """
        Set import mode for headers (data labels)
        """
        text = self.comboBoxHeaders.currentText()
        self.mdl.setTableProperty(self.table_id, 'headers_def', str(text))

        self.lineEditHeadersList.setReadOnly(text == 'file')

        self._update_page_view()

    @Slot(str)
    def slotHeadToImport(self, text):
        """
        Set headers value
        """
        self.mdl.setTableProperty(self.table_id, 'headers_list', str(text))


    @Slot()
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


    @Slot(int)
    def slotHeadLine(self, val):
        """
        Define line from which to import headers
        """
        self.mdl.setTableProperty(self.table_id, 'headers_line', int(val))


    @Slot(int)
    def slotTimeOffsetMode(self, idx):
        """
        Set Offset mode if needed.
        """
        mode = self.comboBoxTimeOffset.currentText()
        self.mdl.setTableProperty(self.table_id, 'time_offset_choice', str(mode))
        if mode == "no":
            self.lineEditTimeOffset.setText('0')
        self._update_page_view()


    @Slot(str)
    def slotTimeOffsetVal(self, val):
        """
        Set Offset mode if needed.
        """
        self.mdl.setTableProperty(self.table_id, 'time_offset_value', val)


    def selectTableFile(self):
        """
        Open a File dialog in order to select a time table.
        """
        _files = []

        title = self.tr("Select input file for time table.")

        #default = os.path.split(self.case['case_path'])[0]
        default = os.path.join(self.case['case_path'], "DATA")

        if hasattr(QFileDialog, 'ReadOnly'):
            options = QFileDialog.Option.DontUseNativeDialog | QFileDialog.Option.ReadOnly
        else:
            options = QFileDialog.Option.DontUseNativeDialog

        filetypes = self.tr('csv file (*.csv);;text file (*.txt);;data file (*.dat);;All files (*)')

        dialog = QFileDialog()
        dialog.setWindowTitle(title)
        dialog.setDirectory(default)
        dialog.setNameFilter(filetypes)

        if hasattr(dialog, 'setOptions'):
            dialog.setOptions(options)
        dialog.setFileMode(QFileDialog.FileMode.ExistingFiles)

        if dialog.exec() == 1:
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


    @Slot()
    def slotPlotTable(self):
        """
        Plot selected time table variables as a function of time.
        """
        if self.table_id < 0:
            return

        fname = self.mdl.getTableFileName(self.table_id)
        table_name = self.mdl.getTableName(self.table_id)

        # Resolve file path
        filepath = fname
        if not os.path.isabs(filepath):
            datap = os.path.join(self.case['case_path'], 'DATA', fname)
            if os.path.isfile(datap):
                filepath = datap

        if not os.path.isfile(filepath):
            QMessageBox.warning(
                self,
                self.tr("File Not Found"),
                self.tr("Time table file '%s' was not found.") % fname
            )
            return

        delimiter = self.mdl.getTableDelimiter(self.table_id)
        delim = delimiter if (delimiter and len(delimiter.strip()) > 0) else None
        if delim and delim.strip():
            delim = delim.strip()

        try:
            skip_rows = int(self._get_mdl_property('skip_rows'))
        except (ValueError, TypeError):
            skip_rows = 0

        cols2import = self._get_mdl_property('cols2import')
        col_ids_str = self._get_mdl_property('col_ids')

        col_indices = None
        if cols2import != 'all' and col_ids_str:
            try:
                # 1-based indexing in GUI
                col_indices = [int(c.strip()) - 1 for c in col_ids_str.split(',') if c.strip()]
            except ValueError:
                col_indices = None

        headers_def = self._get_mdl_property('headers_def')
        headers_list_str = self._get_mdl_property('headers_list')

        labels = []
        if headers_def == 'user' and headers_list_str:
            labels = [h.strip() for h in headers_list_str.split(',') if h.strip()]
        elif headers_def == 'import':
            hl = self._get_headers_from_file(filepath, delimiter)
            if hl:
                labels = [h.strip() for h in hl.split(',') if h.strip()]
        elif headers_list_str:
            labels = [h.strip() for h in headers_list_str.split(',') if h.strip()]

        # Read and parse numerical lines
        try:
            with open(filepath, 'r') as fp:
                all_lines = fp.readlines()
        except Exception as e:
            QMessageBox.warning(
                self,
                self.tr("Error Reading File"),
                self.tr("Could not read file '%s':\n%s") % (fname, str(e))
            )
            return

        data_rows = []
        start_idx = max(0, skip_rows)
        for line in all_lines[start_idx:]:
            sline = line.strip()
            if not sline or sline.startswith('#') or sline.startswith('!'):
                continue

            if delim:
                tokens = [t.strip() for t in sline.split(delim)]
            else:
                tokens = sline.split()

            row_floats = []
            valid_row = True
            for tok in tokens:
                try:
                    row_floats.append(float(tok))
                except ValueError:
                    valid_row = False
                    break

            if not valid_row or not row_floats:
                continue

            if col_indices:
                filtered_row = []
                for idx in col_indices:
                    if 0 <= idx < len(row_floats):
                        filtered_row.append(row_floats[idx])
                if filtered_row:
                    data_rows.append(filtered_row)
            else:
                data_rows.append(row_floats)

        if not data_rows:
            QMessageBox.warning(
                self,
                self.tr("No Data"),
                self.tr("No numeric data rows found in '%s'.") % fname
            )
            return

        num_cols = len(data_rows[0])
        if num_cols < 2:
            QMessageBox.warning(
                self,
                self.tr("Insufficient Columns"),
                self.tr("At least 2 columns (time and at least one variable) are required to plot.")
            )
            return

        # Separate time and variables
        times = [row[0] for row in data_rows]

        time_offset_choice = self._get_mdl_property('time_offset_choice')
        try:
            time_offset_value = float(self._get_mdl_property('time_offset_value'))
        except (ValueError, TypeError):
            time_offset_value = 0.0

        if time_offset_choice == 'yes':
            times = [t - time_offset_value for t in times]

        # Determine variable labels
        time_label = "Time (s)"
        if labels and len(labels) == num_cols:
            time_label = labels[0]
            var_names = labels[1:]
        elif labels and len(labels) == num_cols - 1:
            var_names = labels
        else:
            var_names = [labels[i] if i < len(labels) else f"Variable {i}" for i in range(1, num_cols)]

        variables_dict = {}
        for var_idx, var_name in enumerate(var_names):
            col_idx = var_idx + 1
            if col_idx < num_cols:
                variables_dict[var_name] = [row[col_idx] for row in data_rows]

        dialog = TimeTablePlotDialog(
            parent=self,
            table_name=table_name,
            times=times,
            variables_dict=variables_dict,
            time_label=time_label
        )
        dialog.exec()
