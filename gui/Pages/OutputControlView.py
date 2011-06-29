# -*- coding: utf-8 -*-
#
#-------------------------------------------------------------------------------
#
#     This file is part of the Code_Saturne User Interface, element of the
#     Code_Saturne CFD tool.
#
#     Copyright (C) 1998-2011 EDF S.A., France
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
This module manages the layout of outputs control:
- listing printing
- post-processing and relationship with the FVM library
- monitoring points

This module defines the following classes:
- StandardItemModelMonitoring
- MonitoringPointDelegate
- OutputControliew
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import string
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
from Base.Common import LABEL_LENGTH_MAX
from OutputControlForm import Ui_OutputControlForm
import Base.QtPage as QtPage
from Base.QtPage import ComboModel, DoubleValidator, RegExpValidator, setGreenColor
from Pages.OutputControlModel import OutputControlModel
from Pages.ConjugateHeatTransferModel import ConjugateHeatTransferModel
from Pages.MobileMeshModel import MobileMeshModel
from Pages.QMeiEditorView import QMeiEditorView

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("OutputControlView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Line edit delegate for the label Writer
#-------------------------------------------------------------------------------

class LabelWriterDelegate(QItemDelegate):
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

            if new_plabel in model.mdl.getWriterLabelList():
                default = {}
                default['label']  = self.old_plabel
                default['list']   = model.mdl.getWriterLabelList()
                default['regexp'] = self.regExp
                log.debug("setModelData -> default = %s" % default)

                from Pages.VerifyExistenceLabelDialogView import VerifyExistenceLabelDialogView
                dialog = VerifyExistenceLabelDialogView(self.parent, default)
                if dialog.exec_():
                    result = dialog.get_result()
                    new_plabel = result['label']
                    log.debug("setModelData -> result = %s" % result)
                else:
                    new_plabel = self.old_plabel

            model.setData(index, QVariant(QString(new_plabel)), Qt.DisplayRole)

#-------------------------------------------------------------------------------
# Line edit delegate for the label Mesh
#-------------------------------------------------------------------------------

class LabelMeshDelegate(QItemDelegate):
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

            if new_plabel in model.mdl.getMeshLabelList():
                default = {}
                default['label']  = self.old_plabel
                default['list']   = model.mdl.getMeshLabelList()
                default['regexp'] = self.regExp
                log.debug("setModelData -> default = %s" % default)

                from Pages.VerifyExistenceLabelDialogView import VerifyExistenceLabelDialogView
                dialog = VerifyExistenceLabelDialogView(self.parent, default)
                if dialog.exec_():
                    result = dialog.get_result()
                    new_plabel = result['label']
                    log.debug("setModelData -> result = %s" % result)
                else:
                    new_plabel = self.old_plabel

            model.setData(index, QVariant(QString(new_plabel)), Qt.DisplayRole)

#-------------------------------------------------------------------------------
# Combo box delegate for the writer format
#-------------------------------------------------------------------------------

class FormatWriterDelegate(QItemDelegate):
    """
    Use of a combo box in the table.
    """
    def __init__(self, parent=None, xml_model=None):
        super(FormatWriterDelegate, self).__init__(parent)
        self.parent = parent
        self.mdl = xml_model #à revoir


    def createEditor(self, parent, option, index):
        editor = QComboBox(parent)
        editor.addItem(QString("EnSight"))
        editor.addItem(QString("MED"))
        editor.addItem(QString("CGNS"))
        editor.installEventFilter(self)
        return editor


    def setEditorData(self, comboBox, index):
        dico = {"ensight": 0, "med": 1, "cgns": 2}
        row = index.row()
        string = index.model().dataWriter[row]['format']
        idx = dico[string]
        comboBox.setCurrentIndex(idx)


    def setModelData(self, comboBox, model, index):
        value = comboBox.currentText()
        selectionModel = self.parent.selectionModel()
        for idx in selectionModel.selectedIndexes():
            if idx.column() == index.column():
                model.setData(idx, QVariant(value))

#-------------------------------------------------------------------------------
# Combo box delegate for the writer format
#-------------------------------------------------------------------------------

class TypeMeshDelegate(QItemDelegate):
    """
    Use of a combo box in the table.
    """
    def __init__(self, parent=None, xml_model=None):
        super(TypeMeshDelegate, self).__init__(parent)
        self.parent = parent
        self.mdl = xml_model #à revoir


    def createEditor(self, parent, option, index):
        editor = QComboBox(parent)
        editor.addItem(QString("cells"))
        editor.addItem(QString("interior faces"))
        editor.addItem(QString("boundary faces"))
        editor.installEventFilter(self)
        return editor


    def setEditorData(self, comboBox, index):
        dico = {"cells": 0, "interior_faces": 1, "boundary_faces": 2}
        row = index.row()
        string = index.model().dataMesh[row]['type']
        idx = dico[string]
        comboBox.setCurrentIndex(idx)


    def setModelData(self, comboBox, model, index):
        value = comboBox.currentText()
        selectionModel = self.parent.selectionModel()
        for idx in selectionModel.selectedIndexes():
            if idx.column() == index.column():
                model.setData(idx, QVariant(value))

#-------------------------------------------------------------------------------
# QLineEdit delegate for localization
#-------------------------------------------------------------------------------

class LocationSelectorDelegate(QItemDelegate):
    def __init__(self, parent, mdl):
        super(LocationSelectorDelegate, self).__init__(parent)
        self.parent = parent
        self.mdl = mdl


    def createEditor(self, parent, option, index):
        editor = QLineEdit(parent)
        return editor


    def setEditorData(self, editor, index):
        self.value = index.model().data(index, Qt.DisplayRole).toString()
        editor.setText(self.value)


    def setModelData(self, editor, model, index):
        value = editor.text()

        if str(value) == "" :
           title = self.tr("Warning")
           msg   = self.tr("Please give a localization")
           QMessageBox.information(self.parent, title, msg)
           return

        if str(value) != "" :
            model.setData(index, QVariant(value), Qt.DisplayRole)



#-------------------------------------------------------------------------------
# Combo box delegate for the variance
#-------------------------------------------------------------------------------

class AssociatedWriterDelegate(QItemDelegate):
    """
    Use of a combo box in the table.
    """
    def __init__(self, parent):
        super(AssociatedWriterDelegate, self).__init__(parent)
        self.parent   = parent


    def createEditor(self, parent, option, index):
        editor = QComboBox(parent)
        self.modelCombo = ComboModel(editor, 1, 1)
        editor.installEventFilter(self)
        return editor


    def setEditorData(self, editor, index):
        l1 = index.model().mdl.getWriterLabelList()
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
# QStandardItemModel for Mesh QTableView
#-------------------------------------------------------------------------------

class StandardItemModelMesh(QStandardItemModel):
    def __init__(self, mdl):
        """
        """
        QStandardItemModel.__init__(self)

        self.setColumnCount(4)
        self.dataMesh = []
        self.mdl = mdl
        self.defaultItem = []
        self.populateModel()

    def populateModel(self):
        self.dicoV2M= {"cells": 'cells',"interior faces" : 'interior_faces', "boundary faces": 'boundary_faces'}
        self.dicoM2V= {"cells" : 'cells',"interior_faces" : 'interior faces', "boundary_faces": 'boundary faces'}
        for id in self.mdl.getMeshIdList():
            row = self.rowCount()
            self.setRowCount(row + 1)

            dico           = {}
            dico['name']  = self.mdl.getMeshLabel(id)
            dico['id'] = id
            dico['type'] = self.mdl.getMeshType(id)
            dico['location'] = self.mdl.getMeshLocation(id)

            self.dataMesh.append(dico)
            if int(id) < 0:
                self.defaultItem.append(row)
            log.debug("populateModel-> dataSolver = %s" % dico)


    def data(self, index, role):
        if not index.isValid():
            return QVariant()

        if role == Qt.DisplayRole:

            row = index.row()
            col = index.column()
            dico = self.dataMesh[row]

            if index.column() == 0:
                return QVariant(dico['name'])
            elif index.column() == 1:
                return QVariant(dico['id'])
            elif index.column() == 2:
                return QVariant(self.dicoM2V[dico['type']])
            elif index.column() == 3:
                return QVariant(dico['location'])
            else:
                return QVariant()

        elif role == Qt.TextAlignmentRole:
            return QVariant(Qt.AlignCenter)

        return QVariant()


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled
        # default item
        if index.row() in self.defaultItem:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable
        if index.column() == 1:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable
        else:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable


    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            if section == 0:
                return QVariant(self.tr("Name"))
            elif section == 1:
                return QVariant(self.tr("Id"))
            elif section == 2:
                return QVariant(self.tr("Type"))
            elif section == 3:
                return QVariant(self.tr("Selection"))
        return QVariant()


    def setData(self, index, value, role=None):


        # Update the row in the table
        row = index.row()
        col = index.column()

        # Label
        if col == 0:
            old_plabel = self.dataMesh[row]['name']
            new_plabel = str(value.toString())
            self.dataMesh[row][col] = new_plabel
            self.mdl.setMeshLabel(str(self.dataMesh[row]['id']), new_plabel)

        if index.column() == 2:
            self.dataMesh[row]['type'] = self.dicoV2M[str(value.toString())]
            self.mdl.setMeshType(self.dataMesh[row]['id'], self.dataMesh[row]['type'])

        if index.column() == 3:
            new_nature = str(value.toString())
            self.dataMesh[row]['location'] = new_nature

        self.emit(SIGNAL("dataChanged(const QModelIndex &, const QModelIndex &)"), index, index)
        return True


    def newData(self, name, mesh_id, mesh_type, location):
        """
        Add a new 'item' into the table.
        """
        dico = {}
        dico['name'] = name
        dico['id'] = mesh_id
        dico['type'] = mesh_type
        dico['location'] = location
        self.dataMesh.append(dico)

        row = self.rowCount()
        self.setRowCount(row + 1)


    def getItem(self, row):
        return self.dataMesh[row]


    def getData(self, row, column):
        return self.dataMesh[row][column]


    def deleteAllData(self):
        """
        Destroy the contents of the list.
        """
        self.dataMesh = []
        self.setRowCount(0)

#-------------------------------------------------------------------------------
# QStandardItemModel for Mesh QTableView
#-------------------------------------------------------------------------------

class StandardItemModelWriter(QStandardItemModel):
    def __init__(self, mdl):
        """
        """
        QStandardItemModel.__init__(self)

        self.setColumnCount(4)
        self.dataWriter = []
        self.mdl = mdl
        self.defaultItem = []
        self.populateModel()

    def populateModel(self):
        self.dicoV2M= {"EnSight": 'ensight',"MED" : 'med', "CGNS": 'cgns'}
        self.dicoM2V= {"ensight" : 'EnSight',"med" : 'MED', "cgns": 'CGNS'}
        for id in self.mdl.getWriterIdList():
            row = self.rowCount()
            self.setRowCount(row + 1)

            dico           = {}
            dico['name']  = self.mdl.getWriterLabel(id)
            dico['id'] = id
            dico['format'] = self.mdl.getWriterFormat(id)
            dico['repertory']  = self.mdl.getWriterRepertory(id)

            self.dataWriter.append(dico)
            if int(id)<0:
                self.defaultItem.append(row)
            log.debug("populateModel-> dataSolver = %s" % dico)


    def data(self, index, role):
        if not index.isValid():
            return QVariant()

        if role == Qt.DisplayRole:

            row = index.row()
            col = index.column()
            dico = self.dataWriter[row]

            if index.column() == 0:
                return QVariant(dico['name'])
            elif index.column() == 1:
                return QVariant(dico['id'])
            elif index.column() == 2:
                return QVariant(self.dicoM2V[dico['format']])
            elif index.column() == 3:
                return QVariant(dico['repertory'])
            else:
                return QVariant()

        elif role == Qt.TextAlignmentRole:
            return QVariant(Qt.AlignCenter)

        return QVariant()


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled
        # default item
        if index.row() in self.defaultItem:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable
        if index.column() == 1:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable
        else:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable


    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            if section == 0:
                return QVariant(self.tr("Name"))
            elif section == 1:
                return QVariant(self.tr("Id"))
            elif section == 2:
                return QVariant(self.tr("Format"))
            elif section == 1:
                return QVariant(self.tr("Repertory"))
        return QVariant()


    def setData(self, index, value, role=None):


        # Update the row in the table
        row = index.row()
        col = index.column()

        # Label
        if col == 0:
            old_plabel = self.dataWriter[row]['name']
            new_plabel = str(value.toString())
            self.dataWriter[row]['name'] = new_plabel
            self.mdl.setWriterLabel(self.dataWriter[row]['id'], new_plabel)

        if index.column() == 2:
            self.dataWriter[row]['format'] = self.dicoV2M[str(value.toString())]
            self.mdl.setWriterFormat(self.dataWriter[row]['id'], self.dataWriter[row]['format'])

        if col == 3:
            old_rep = self.dataWriter[row]['repertory']
            new_rep = str(value.toString())
            self.dataWriter[row]['repertory'] = new_rep
            self.mdl.setWriterRepertory(self.dataWriter[row]['id'], new_rep)


        self.emit(SIGNAL("dataChanged(const QModelIndex &, const QModelIndex &)"), index, index)
        return True


    def newData(self, name, writer_id, writer_format, repertory):
        """
        Add a new 'item' into the table.
        """
        dico = {}
        dico['name'] = name
        dico['id'] = writer_id
        dico['format'] = writer_format
        dico['repertory'] = repertory
        self.dataWriter.append(dico)

        row = self.rowCount()
        self.setRowCount(row + 1)


    def getItem(self, row):
        return self.dataWriter[row]


    def getData(self, row, column):
        return self.dataWriter[row][column]


    def deleteAllData(self):
        """
        Destroy the contents of the list.
        """
        self.dataWriter = []
        self.setRowCount(0)



#-------------------------------------------------------------------------------
# StandarItemModel class
#-------------------------------------------------------------------------------

class StandardItemModelAssociatedWriter(QStandardItemModel):
    """
    """
    def __init__(self, parent, mdl, mesh_id):
        """
        """
        QStandardItemModel.__init__(self)

        self.headers = [self.tr("Writer")]

        self.setColumnCount(len(self.headers))

        self._data = []
        self.parent = parent
        self.mdl  = mdl
        self.mesh_id = mesh_id


    def data(self, index, role):
        if not index.isValid():
            return QVariant()

        row = index.row()
        col = index.column()

        if role == Qt.DisplayRole:
            return QVariant(self._data[row])

        return QVariant()


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled
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
            writer = str(value.toString())
            self._data[row] = writer
            writer_list = []
            for r in range(self.rowCount()):
                writer_list.append(self.mdl.getWriterIdFromLabel(self._data[r]))
            self.mdl.setAssociatedWriterChoice(self.mesh_id, writer_list)

        self.emit(SIGNAL("dataChanged(const QModelIndex &, const QModelIndex &)"), index, index)
        return True


    def getData(self, index):
        row = index.row()
        return self._data[row]


    def newItem(self, label):
        """
        Add an item in the table view
        """
        if not self.mdl.getWriterIdList():
            title = self.tr("Warning")
            msg   = self.tr("There is no writer.\n"\
                            "Please define a writer.")
            QMessageBox.warning(self.parent, title, msg)
            return
        row = self.rowCount()
        self.setRowCount(row+1)
        self._data.append(label)


    def getItem(self, row):
        """
        Return the values for an item.
        """
        var = self._data[row]
        return var


    def deleteItem(self, row):
        """
        Delete the row in the model.
        """
        log.debug("deleteItem row = %i " % row)

        del self._data[row]
        row = self.rowCount()
        self.setRowCount(row-1)


    def deleteAllData(self):
        """
        Destroy the contents of the list.
        """
        self._data = []
        self.setRowCount(0)





#-------------------------------------------------------------------------------
# QStandardItemModel for monitoring points QTableView
#-------------------------------------------------------------------------------

class StandardItemModelMonitoring(QStandardItemModel):
    def __init__(self):
        """
        """
        QStandardItemModel.__init__(self)

        self.setColumnCount(4)
        self.dataMonitoring = []


    def data(self, index, role):
        if not index.isValid():
            return QVariant()

        if role == Qt.DisplayRole:

            row = index.row()
            dico = self.dataMonitoring[row]

            if index.column() == 0:
                return QVariant(dico['n'])
            elif index.column() == 1:
                return QVariant(dico['X'])
            elif index.column() == 2:
                return QVariant(dico['Y'])
            elif index.column() == 3:
                return QVariant(dico['Z'])
            else:
                return QVariant()

        elif role == Qt.TextAlignmentRole:
            return QVariant(Qt.AlignCenter)

        return QVariant()


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled
        if index.column() == 0:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable
        else:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable


    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            if section == 0:
                return QVariant(self.tr("n"))
            elif section == 1:
                return QVariant(self.tr("X"))
            elif section == 2:
                return QVariant(self.tr("Y"))
            elif section == 3:
                return QVariant(self.tr("Z"))
        return QVariant()


    def setData(self, index, value, role=None):
        row = index.row()
        if index.column() == 0:
            n, ok = value.toInt()
            self.dataMonitoring[row]['n'] = n
        elif index.column() == 1:
            X, ok = value.toDouble()
            self.dataMonitoring[row]['X'] = X
        elif index.column() == 2:
            Y, ok = value.toDouble()
            self.dataMonitoring[row]['Y'] = Y
        elif index.column() == 3:
            Z, ok = value.toDouble()
            self.dataMonitoring[row]['Z'] = Z

        self.emit(SIGNAL("dataChanged(const QModelIndex &, const QModelIndex &)"), index, index)
        return True


    def insertData(self, num, X, Y, Z):
        """
        Add a new 'item' into the table.
        """
        dico = {}
        dico['n'] = num
        dico['X'] = X
        dico['Y'] = Y
        dico['Z'] = Z
        self.dataMonitoring.append(dico)

        row = self.rowCount()
        self.setRowCount(row + 1)


    def deleteAllData(self):
        """
        Destroy the contents of the list.
        """
        self.dataMonitoring = []
        self.setRowCount(0)

#-------------------------------------------------------------------------------
# QItemDelegate for monitoring points QTableView
#-------------------------------------------------------------------------------

class MonitoringPointDelegate(QItemDelegate):
    def __init__(self, parent=None, xml_model=None):
        """ Construtor.

        @param: parent ancestor object
        @xml_model: monitoring points model
        """
        super(MonitoringPointDelegate, self).__init__(parent)
        self.table = parent
        self.mdl = xml_model


    def createEditor(self, parent, option, index):
        if index.column() == 0:
            editor = QFrame(parent)
        else:
            editor = QLineEdit(parent)
            editor.setValidator(QtPage.DoubleValidator(editor))
            editor.setFrame(False)
            self.connect(editor, SIGNAL("returnPressed()"), self.commitAndCloseEditor)
            editor.setCursorPosition(0)
        return editor


    def commitAndCloseEditor(self):
        editor = self.sender()
        if isinstance(editor, QLineEdit):
            self.emit(SIGNAL("commitData(QWidget*)"), editor)
            self.emit(SIGNAL("closeEditor(QWidget*)"), editor)


    def setEditorData(self, editor, index):
        text = index.model().data(index, Qt.DisplayRole).toString()
        if isinstance(editor, QLineEdit):
            editor.setText(text)


    def setModelData(self, editor, model, index):
        if isinstance(editor, QLineEdit):
            if not editor.isModified():
                return

            item = editor.text()
            selectionModel = self.table.selectionModel()
            for index in selectionModel.selectedRows(index.column()):
                model.setData(index, QVariant(item), Qt.DisplayRole)
                dico = model.dataMonitoring[index.row()]
                self.mdl.replaceMonitoringPointCoordinates(str(dico['n']),
                                                           float(dico['X']),
                                                           float(dico['Y']),
                                                           float(dico['Z']))

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class OutputControlView(QWidget, Ui_OutputControlForm):
    """
    """
    def __init__(self, parent, case, tree):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_OutputControlForm.__init__(self)
        self.setupUi(self)

        self.browser = tree
        self.case = case
        self.mdl = OutputControlModel(self.case)

        # Combo models

        self.modelOutput         = QtPage.ComboModel(self.comboBoxOutput,3,1)
        self.modelFrequency = QtPage.ComboModel(self.comboBoxFrequency,5,1)
        self.modelTimeDependency         = QtPage.ComboModel(self.comboBoxTimeDependency,3,1)
        self.modelFormat         = QtPage.ComboModel(self.comboBoxFormat,2,1)
        self.modelPolygon        = QtPage.ComboModel(self.comboBoxPolygon,3,1)
        self.modelPolyhedra      = QtPage.ComboModel(self.comboBoxPolyhedra,3,1)
        self.modelHisto          = QtPage.ComboModel(self.comboBoxHisto,3,1)
        self.modelProbeFmt       = QtPage.ComboModel(self.comboBoxProbeFmt,2,1)

        self.modelOutput.addItem(self.tr("No output"), 'None')
        self.modelOutput.addItem(self.tr("Output listing at each time step"), 'At each step')
        self.modelOutput.addItem(self.tr("Output every 'n' time steps"), 'Frequency_l')

        self.modelFrequency.addItem(self.tr("Only at the end of calculation"), 'end')
        self.modelFrequency.addItem(self.tr("At each time step"), 'time')
        self.modelFrequency.addItem(self.tr("Post-processing every 'n' time steps"), 'time_steps')
        self.modelFrequency.addItem(self.tr("Post-processing every 'x' second(s)"), 'second')
        self.modelFrequency.addItem(self.tr("Using a formula"), 'formula')

        self.modelTimeDependency.addItem(self.tr("Fixed mesh"), 'fixed_mesh')
        self.modelTimeDependency.addItem(self.tr("Transient coordinates"), 'transient_coordinates')
        self.modelTimeDependency.addItem(self.tr("Transient connectivity"), 'transient_connectivity')

        #ale = self.ale()

        self.modelFormat.addItem(self.tr("binary"), 'binary')
        self.modelFormat.addItem(self.tr("text"), 'text')

        self.modelPolygon.addItem(self.tr("display"), 'display')
        self.modelPolygon.addItem(self.tr("discard"), 'discard_polygons')
        self.modelPolygon.addItem(self.tr("subdivide"), 'divide_polygons')

        self.modelPolyhedra.addItem(self.tr("display"), 'display')
        self.modelPolyhedra.addItem(self.tr("discard"), 'discard_polyhedra')
        self.modelPolyhedra.addItem(self.tr("subdivide"), 'divide_polyhedra')

        self.modelHisto.addItem(self.tr("No monitoring points file"), 'None')
        self.modelHisto.addItem(self.tr("Monitoring points files at each time step"), 'At each step')
        self.modelHisto.addItem(self.tr("Monitoring points file every 'n' time steps"), 'Frequency_h')
        self.modelHisto.addItem(self.tr("Monitoring points file every 'x' second(s)"), 'Frequency_h_x')

        self.modelProbeFmt.addItem(self.tr(".dat"), 'DAT')
        self.modelProbeFmt.addItem(self.tr(".csv"), 'CSV')

        # Hide time frequency (in s) when calculation is steady
        if self.isSteady() != 1:
            self.modelHisto.disableItem(3)

        # Model for QTableView

        self.modelMonitoring = StandardItemModelMonitoring()
        self.tableViewPoints.setModel(self.modelMonitoring)
        self.tableViewPoints.resizeColumnToContents(0)
        self.tableViewPoints.resizeRowsToContents()
        self.tableViewPoints.setAlternatingRowColors(True)
        self.tableViewPoints.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.tableViewPoints.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.tableViewPoints.setEditTriggers(QAbstractItemView.DoubleClicked)
        self.tableViewPoints.horizontalHeader().setResizeMode(QHeaderView.Stretch)
        delegate = MonitoringPointDelegate(self.tableViewPoints, self.mdl)
        self.tableViewPoints.setItemDelegate(delegate)

        self.modelWriter = StandardItemModelWriter(self.mdl)
        self.tableViewWriter.setModel(self.modelWriter)
        self.tableViewWriter.resizeColumnToContents(0)
        self.tableViewWriter.resizeRowsToContents()
        self.tableViewWriter.setAlternatingRowColors(True)
        self.tableViewWriter.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.tableViewWriter.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.tableViewWriter.setEditTriggers(QAbstractItemView.DoubleClicked)
        self.tableViewWriter.horizontalHeader().setResizeMode(QHeaderView.Stretch)

        delegate_label_writer = LabelWriterDelegate(self.tableViewWriter)
        self.tableViewWriter.setItemDelegateForColumn(0, delegate_label_writer)
        self.tableViewWriter.setItemDelegateForColumn(3, delegate_label_writer)
        delegate_format = FormatWriterDelegate(self.tableViewWriter, self.mdl)
        self.tableViewWriter.setItemDelegateForColumn(2, delegate_format)

        self.modelMesh = StandardItemModelMesh(self.mdl)
        self.tableViewMesh.setModel(self.modelMesh)
        self.tableViewMesh.resizeColumnToContents(0)
        self.tableViewMesh.resizeRowsToContents()
        self.tableViewMesh.setAlternatingRowColors(True)
        self.tableViewMesh.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.tableViewMesh.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.tableViewMesh.setEditTriggers(QAbstractItemView.DoubleClicked)
        self.tableViewMesh.horizontalHeader().setResizeMode(QHeaderView.Stretch)

        delegate_label_mesh = LabelMeshDelegate(self.tableViewMesh)
        self.tableViewMesh.setItemDelegateForColumn(0, delegate_label_mesh)
        delegate_type = TypeMeshDelegate(self.tableViewMesh, self.mdl)
        self.tableViewMesh.setItemDelegateForColumn(2, delegate_type)
        delegate_location = LocationSelectorDelegate(self.tableViewMesh, self.mdl)
        self.tableViewMesh.setItemDelegateForColumn(3, delegate_location)

        # Connections

        self.connect(self.modelWriter,     SIGNAL("dataChanged(const QModelIndex &, const QModelIndex &)"), self.dataChanged)
        self.connect(self.tableViewMesh, SIGNAL("clicked(const QModelIndex &)"), self.slotSelectMesh)
        self.connect(self.tableViewWriter, SIGNAL("clicked(const QModelIndex &)"), self.slotSelectWriter)
        self.connect(self.comboBoxOutput, SIGNAL("activated(const QString&)"), self.slotOutputListing)
        self.connect(self.comboBoxTimeDependency, SIGNAL("activated(const QString&)"), self.slotWriterTimeDependency)
        self.connect(self.checkBoxAllVariables, SIGNAL("clicked()"), self.slotAllVariables)
        self.connect(self.lineEditNTLIST, SIGNAL("textChanged(const QString &)"), self.slotListingFrequency)
        self.connect(self.comboBoxFrequency, SIGNAL("activated(const QString&)"), self.slotWriterFrequencyChoice)
        self.connect(self.lineEditFrequency, SIGNAL("textChanged(const QString &)"), self.slotWriterFrequency)
        self.connect(self.lineEditFrequencyTime, SIGNAL("textChanged(const QString &)"), self.slotWriterFrequencyTime)
        self.connect(self.comboBoxFormat, SIGNAL("activated(const QString&)"), self.slotWriterOptions)
        self.connect(self.comboBoxPolygon, SIGNAL("activated(const QString&)"), self.slotWriterOptions)
        self.connect(self.comboBoxPolyhedra, SIGNAL("activated(const QString&)"), self.slotWriterOptions)
        self.connect(self.checkBoxBigEndian, SIGNAL("clicked()"), self.slotWriterOptions)
        self.connect(self.pushButtonFrequency, SIGNAL("clicked()"), self.slotWriterFrequencyFormula)

        self.connect(self.pushButtonAddWriter, SIGNAL("clicked()"), self.slotAddWriter)
        self.connect(self.pushButtonDeleteWriter, SIGNAL("clicked()"), self.slotDeleteWriter)
        self.connect(self.pushButtonAddMesh, SIGNAL("clicked()"), self.slotAddMesh)
        self.connect(self.pushButtonDeleteMesh, SIGNAL("clicked()"), self.slotDeleteMesh)
        self.connect(self.pushButtonAddAssociatedWriter, SIGNAL("clicked()"), self.slotAddAssociatedWriter)
        self.connect(self.pushButtonDeleteAssociatedWriter, SIGNAL("clicked()"), self.slotDeleteAssociatedWriter)
        self.connect(self.toolButtonAdd, SIGNAL("clicked()"), self.slotAddMonitoringPoint)
        self.connect(self.toolButtonDelete, SIGNAL("clicked()"), self.slotDeleteMonitoringPoints)
        self.connect(self.comboBoxHisto, SIGNAL("activated(const QString&)"), self.slotMonitoringPoint)
        self.connect(self.lineEditHisto, SIGNAL("textChanged(const QString &)"), self.slotMonitoringPointFrequency)
        self.connect(self.lineEditFRHisto, SIGNAL("textChanged(const QString &)"), self.slotMonitoringPointFrequencyTime)
        self.connect(self.comboBoxProbeFmt, SIGNAL("activated(const QString&)"), self.slotOutputProbeFmt)

        # Validators

        validatorNTLIST = QtPage.IntValidator(self.lineEditNTLIST, min=1)
        validatorFrequency  = QtPage.IntValidator(self.lineEditFrequency, min=1)
        validatorNTHIST = QtPage.IntValidator(self.lineEditHisto, min=1)
        validatorFrequencyTime  = QtPage.DoubleValidator(self.lineEditFrequencyTime)
        validatorFRHIST = QtPage.DoubleValidator(self.lineEditFRHisto)
        self.lineEditNTLIST.setValidator(validatorNTLIST)
        self.lineEditFrequency.setValidator(validatorFrequency)
        self.lineEditHisto.setValidator(validatorNTHIST)
        self.lineEditFrequencyTime.setValidator(validatorFrequencyTime)
        self.lineEditFRHisto.setValidator(validatorFRHIST)

        # Initialisation of the listing frequency

        ntlist = self.mdl.getListingFrequency()
        if ntlist == -1:
            m = "None"
        elif ntlist == 1:
            m = "At each step"
        else:
            m = "Frequency_l"
        self.modelOutput.setItem(str_model=m)
        t = self.modelOutput.dicoM2V[m]
        self.lineEditNTLIST.setText(QString(str(ntlist)))
        self.slotOutputListing(t)


        # Initialisation of the monitoring points files

        m = self.mdl.getMonitoringPointType()
        if m == 'Frequency_h_x' :
            frhist = self.mdl.getMonitoringPointFrequencyTime()
            self.lineEditFRHisto.setText(QString(str(frhist)))
        else :
            nthist = self.mdl.getMonitoringPointFrequency()
            self.lineEditHisto.setText(QString(str(nthist)))
        self.modelHisto.setItem(str_model=m)
        t = self.modelHisto.dicoM2V[m]
        self.slotMonitoringPoint(t)


        # Monitoring points initialisation

        for n in range(self.mdl.getNumberOfMonitoringPoints()):
            name = str(n+1)
            X, Y, Z = self.mdl.getMonitoringPointCoordinates(name)
            self.__insertMonitoringPoint(name, X, Y, Z)

        # Writer initialisation

        self.groupBoxFrequency.hide()
        self.groupBoxTimeDependency.hide()
        self.groupBoxOptions.hide()

        # Mesh initialisation

        self.groupBoxVariable.hide()
        self.groupBoxAssociatedWriter.hide()

        # values of probes format

        fmt = self.mdl.getMonitoringPointFormat()
        self.modelProbeFmt.setItem(str_model=fmt)


    @pyqtSignature("const QString &")
    def slotOutputListing(self, text):
        """
        INPUT choice of the output listing
        """
        listing = self.modelOutput.dicoV2M[str(text)]
        log.debug("slotOutputListing-> listing = %s" % listing)

        if listing == "None":
            ntlist = -1
            self.mdl.setListingFrequency(ntlist)
            self.lineEditNTLIST.setText(QString(str(ntlist)))
            self.lineEditNTLIST.setDisabled(True)

        elif listing == "At each step":
            ntlist = 1
            self.lineEditNTLIST.setText(QString(str(ntlist)))
            self.lineEditNTLIST.setDisabled(True)

        elif listing == "Frequency_l":
            self.lineEditNTLIST.setEnabled(True)
            ntlist, ok = self.lineEditNTLIST.text().toInt()
            if ntlist < 1:
                ntlist = 1
                self.lineEditNTLIST.setText(QString(str(ntlist)))


    @pyqtSignature("const QString &")
    def slotListingFrequency(self, text):
        """
        Input the frequency of the listing output
        """
        n, ok = text.toInt()
        if self.sender().validator().state == QValidator.Acceptable:
            log.debug("slotListingFrequency-> NTLIST = %s" % n)
            self.mdl.setListingFrequency(n)


    def __insertWriter(self, name, writer_id, format, repertory):
        """
        Add a new 'item' into the Hlist.
        """
        self.modelWriter.newData(name, writer_id, format, repertory)


    @pyqtSignature("")
    def slotAddWriter(self):
        """
        Add one monitoring point with these coordinates in the list in the Hlist
        The number of the monitoring point is added at the precedent one
        """
        writer_id = self.mdl.addWriter()
        self.__insertWriter(self.mdl.getWriterLabel(writer_id),
                                     writer_id,
                                     self.mdl.getWriterFormat(writer_id),
                                     self.mdl.getWriterRepertory(writer_id))


    @pyqtSignature("")
    def slotDeleteWriter(self):
        """
        Just delete the current selected entries from the Hlist and
        of course from the XML file.
        """
        list = []
        selectionModel = self.tableViewWriter.selectionModel()
        for index in selectionModel.selectedRows():
            w = self.modelWriter.getItem(index.row())['id']
            if int(w) < 0:
                title = self.tr("Warning")
                msg   = self.tr("You can't delete a default writer")
                QMessageBox.information(self, title, msg)
                return
            list.append(str(w))

        self.mdl.deleteWriter(list)

        self.modelWriter.deleteAllData()
        list_writer = []
        for writer in self.mdl.getWriterIdList():
            if int(writer)>0:
                list_writer.append(writer)
        for writer in self.mdl.getWriterIdList():
            if int(writer)<0:
                label = self.mdl.getWriterLabel(writer)
                format = self.mdl.getWriterFormat(writer)
                repertory = self.mdl.getWriterRepertory(writer)
                self.__insertWriter(label, writer, format, repertory)
        new_id = 0
        for writer in list_writer:
            new_id = new_id + 1
            label = self.mdl.getWriterLabel(writer)
            format = self.mdl.getWriterFormat(writer)
            repertory = self.mdl.getWriterRepertory(writer)
            self.__insertWriter(label, str(new_id), format, repertory)


    @pyqtSignature("const QModelIndex &, const QModelIndex &")
    def dataChanged(self, topLeft, bottomRight):
        for row in range(topLeft.row(), bottomRight.row()+1):
            self.tableViewWriter.resizeRowToContents(row)
        for col in range(topLeft.column(), bottomRight.column()+1):
            self.tableViewWriter.resizeColumnToContents(col)
        cindex = self.tableViewWriter.currentIndex()
        if cindex != (-1,-1):
            row_writer = cindex.row()
            writer_id = self.modelWriter.getItem(row_writer)['id']
            options = self.mdl.getWriterOptions(writer_id)
            self.__updateOptionsFormat(options, row_writer)
            self.showAssociatedWriterTable()
        #print 'passage dans data changed'


    def showAssociatedWriterTable(self):
        cindex = self.tableViewMesh.currentIndex()
        if cindex != (-1,-1):
            row = cindex.row()
            mesh_id = self.modelMesh.getItem(row)['id']

            self.modelAssociatedWriter = StandardItemModelAssociatedWriter(self, self.mdl, mesh_id)
            self.tableViewAssociatedWriter.horizontalHeader().setResizeMode(QHeaderView.Stretch)

            delegate_associated_writer = AssociatedWriterDelegate(self.tableViewAssociatedWriter)
            self.tableViewAssociatedWriter.setItemDelegateForColumn(0, delegate_associated_writer)

            self.tableViewAssociatedWriter.reset()
            self.modelAssociatedWriter = StandardItemModelAssociatedWriter(self, self.mdl, mesh_id)
            self.tableViewAssociatedWriter.setModel(self.modelAssociatedWriter)
            self.modelAssociatedWriter.deleteAllData()
            writer_row = 0
            for n in self.mdl.getAssociatedWriterIdList(mesh_id):
                label = self.mdl.getWriterLabel(n)
                self.__insertAssociatedWriter(label)
                writer_row = writer_row +1


    @pyqtSignature("const QModelIndex&")
    def slotSelectWriter(self, index):
        #model = HeadLossesModel(self.case)
        cindex = self.tableViewWriter.currentIndex()
        if cindex != (-1,-1):
            row = cindex.row()
            writer_id = self.modelWriter.getItem(row)['id']
            self.groupBoxFrequency.show()
            self.groupBoxTimeDependency.show()
            self.groupBoxOptions.show()

            frequency_choice = self.mdl.getWriterFrequencyChoice(writer_id)
            self.modelFrequency.setItem(str_model=frequency_choice)



            if frequency_choice == "end":
                ntchr = -1
                self.mdl.setWriterFrequency(writer_id, ntchr)
                self.lineEditFrequency.hide()
                self.lineEditFrequencyTime.hide()
                self.pushButtonFrequency.setEnabled(False)
                setGreenColor(self.pushButtonFrequency, False)

            if frequency_choice == "time":
                ntchr = 1
                self.mdl.setWriterFrequency(writer_id, ntchr)
                self.lineEditFrequency.hide()
                self.lineEditFrequencyTime.hide()
                self.pushButtonFrequency.setEnabled(False)
                setGreenColor(self.pushButtonFrequency, False)

            if frequency_choice == "time_steps":
                self.lineEditFrequency.show()
                self.lineEditFrequency.setEnabled(True)
                ntchr = self.mdl.getWriterFrequency(writer_id)
                if ntchr < 1:
                    ntchr = 1
                    self.mdl.setWriterFrequency(writer_id, ntchr)
                self.lineEditFrequency.setText(QString(str(ntchr)))
                self.lineEditFrequencyTime.hide()
                self.pushButtonFrequency.setEnabled(False)
                setGreenColor(self.pushButtonFrequency, False)

            if frequency_choice == "second":
                self.lineEditFrequency.hide()
                self.lineEditFrequencyTime.show()
                frchr = self.mdl.getWriterFrequencyTime(writer_id)
                self.lineEditFrequencyTime.setText(QString(str(frchr)))
                self.pushButtonFrequency.setEnabled(False)
                setGreenColor(self.pushButtonFrequency, False)

            if frequency_choice == "formula":
                self.lineEditFrequency.hide()
                self.lineEditFrequencyTime.hide()
                self.pushButtonFrequency.setEnabled(True)
                setGreenColor(self.pushButtonFrequency, True)


            time_dependency = self.mdl.getWriterTimeDependency(writer_id)
            self.modelTimeDependency.setItem(str_model=time_dependency)
            options = self.mdl.getWriterOptions(writer_id)
            self.__updateOptionsFormat(options, row)


    @pyqtSignature("const QString &")
    def slotWriterFrequencyChoice(self, text):
        """
        INPUT choice of the output frequency for a writer
        """
        cindex = self.tableViewWriter.currentIndex()
        if cindex != (-1,-1):
            row = cindex.row()
            writer_id = self.modelWriter.getItem(row)['id']
            chrono = self.modelFrequency.dicoV2M[str(text)]
            log.debug("slotOutputPostpro-> chrono = %s" % chrono)
            self.mdl.setWriterFrequencyChoice(writer_id, chrono)
            if chrono == "end":
                ntchr = -1
                self.mdl.setWriterFrequency(writer_id, ntchr)
                self.lineEditFrequency.hide()
                self.lineEditFrequencyTime.hide()
                self.pushButtonFrequency.setEnabled(False)
                setGreenColor(self.pushButtonFrequency, False)

            elif chrono == "time":
                ntchr = 1
                self.mdl.setWriterFrequency(writer_id, ntchr)
                self.lineEditFrequency.hide()
                self.lineEditFrequencyTime.hide()
                self.pushButtonFrequency.setEnabled(False)
                setGreenColor(self.pushButtonFrequency, False)

            elif chrono == "time_steps":
                self.lineEditFrequency.show()
                self.lineEditFrequency.setEnabled(True)
                self.pushButtonFrequency.setEnabled(False)
                ntchr = self.mdl.getWriterFrequency(writer_id)
                if ntchr < 1:
                    ntchr = 1
                    self.mdl.setWriterFrequency(writer_id, ntchr)
                self.lineEditFrequency.setText(QString(str(ntchr)))
                self.lineEditFrequencyTime.hide()
                setGreenColor(self.pushButtonFrequency, False)

            elif chrono == "second":
                self.lineEditFrequency.hide()
                self.lineEditFrequencyTime.show()
                self.pushButtonFrequency.setEnabled(False)
                frchr = self.mdl.getWriterFrequencyTime(writer_id)
                self.lineEditFrequencyTime.setText(QString(str(frchr)))
                setGreenColor(self.pushButtonFrequency, False)

            elif chrono == "formula":
                self.lineEditFrequency.hide()
                self.lineEditFrequencyTime.hide()
                self.pushButtonFrequency.setEnabled(True)
                setGreenColor(self.pushButtonFrequency, True)


    @pyqtSignature("const QString &")
    def slotWriterFrequency(self, text):
        """
        Input the frequency of the post-processing output
        """
        cindex = self.tableViewWriter.currentIndex()
        if cindex != (-1,-1):
            row = cindex.row()
            writer_id = self.modelWriter.getItem(row)['id']
            self.lineEditFrequency.setEnabled(True)
            n, ok = self.lineEditFrequency.text().toInt()
            if self.sender().validator().state == QValidator.Acceptable:
                log.debug("slotPostproFrequency-> NTCHR = %s" % n)
                self.mdl.setWriterFrequency(writer_id, n)


    @pyqtSignature("const QString &")
    def slotWriterFrequencyTime(self, text):
        """
        Input the frequency of the post-processing output
        """
        cindex = self.tableViewWriter.currentIndex()
        if cindex != (-1,-1):
            row = cindex.row()
            writer_id = self.modelWriter.getItem(row)['id']
            n, ok = text.toDouble()
            if self.sender().validator().state == QValidator.Acceptable:
                log.debug("slotPostproFrequencyTime-> FRCHR = %s" % n)
                self.mdl.setWriterFrequencyTime(writer_id, n)


    @pyqtSignature("")
    def slotWriterFrequencyFormula(self):
        """
        """
        cindex = self.tableViewWriter.currentIndex()
        if cindex != (-1,-1):
            row = cindex.row()
            writer_id = self.modelWriter.getItem(row)['id']
            exp = self.mdl.getWriterFrequencyFormula(writer_id)
            if not exp:
                exp = """iactive = 1;\n"""
            exa = """#example:"""
            req = [('iactive', 'at a time step the writer is active or not')]
            sym = [('t', 'current time'),
                ('niter', 'current time step')]
            dialog = QMeiEditorView(self,expression = exp,
                                        required   = req,
                                        symbols    = sym,
                                        examples   = exa)
            if dialog.exec_():
                result = dialog.get_result()
                log.debug("slotWriterFrequencyFormula -> %s" % str(result))
                self.mdl.setWriterFrequencyFormula(writer_id, result)
                setGreenColor(self.pushButtonFrequency, False)


    @pyqtSignature("const QString &")
    def slotWriterTimeDependency(self, text):
        """
        Input type of post-processing for mesh
        """
        cindex = self.tableViewWriter.currentIndex()
        if cindex != (-1,-1):
            row = cindex.row()
            writer_id = self.modelWriter.getItem(row)['id']
            self.mdl.setWriterTimeDependency(writer_id, self.modelTimeDependency.dicoV2M[str(text)])


    @pyqtSignature("")
    def slotWriterOptions(self):
        """
        Create characters ligne for command of format's options
        """
        cindex = self.tableViewWriter.currentIndex()
        if cindex != (-1,-1):
            row = cindex.row()
            writer_id = self.modelWriter.getItem(row)['id']
            line = []
            opt_format = self.modelFormat.dicoV2M[str(self.comboBoxFormat.currentText())]
            line.append(opt_format)

            if self.checkBoxBigEndian.isChecked():
                line.append('big_endian')

            opt_polygon = self.modelPolygon.dicoV2M[str(self.comboBoxPolygon.currentText())]
            opt_polyhed = self.modelPolyhedra.dicoV2M[str(self.comboBoxPolyhedra.currentText())]
            if opt_polygon != 'display': line.append(opt_polygon)
            if opt_polyhed != 'display': line.append(opt_polyhed)

            l = string.join(line, ',')
            log.debug("slotOutputOptions-> OPTCHR = %s" % l)
            self.mdl.setWriterOptions(writer_id, l)


    def __updateOptionsFormat(self, line, row):# appelé en initialisation
        """
        Update ligne for command of format's options at each modification of
        post processing format
        """
        list = string.split(line, ',')
        format = self.modelWriter.getItem(row)['format']

        # update widgets from the options list

        for opt in list:

            if opt == 'binary' or opt == 'text' :
                self.modelFormat.setItem(str_model=opt)

            if opt == 'discard polygons' or opt == 'divide_polygons':
                self.modelPolygon.setItem(str_model=opt)

            if opt == 'discard polyhedra' or opt == 'divide_polyhedra':
                self.modelPolyhedra.setItem(str_model=opt)

            if format == 'ensight':
                if opt == 'big_endian':
                    self.checkBoxBigEndian.setChecked(True)

        if 'discard_polygons' not in list and 'divide_polygons' not in list:
            self.modelPolygon.setItem(str_model="display")
        if 'discard_polyhedra' not in list and 'divide_polyhedra' not in list:
            self.modelPolyhedra.setItem(str_model="display")
        if 'big_endian' not in list:
            self.checkBoxBigEndian.setChecked(False)

        # enable and disable options related to the format
        if format != "ensight":

            if format == "cgns":
                self.modelPolyhedra.setItem(str_model='divide_polyhedra')
                self.modelPolyhedra.disableItem(str_model='display')

            self.modelFormat.setItem(str_model="binary")
            self.modelFormat.disableItem(str_model='text')
            self.labelBigEndian.setEnabled(False)
            self.checkBoxBigEndian.setEnabled(False)
        else:
            self.modelFormat.enableItem(str_model='text')
            self.comboBoxFormat.setEnabled(True)
            self.labelBigEndian.setEnabled(True)
            self.checkBoxBigEndian.setEnabled(True)
            self.modelPolyhedra.enableItem(str_model='display')
            self.comboBoxPolyhedra.setEnabled(True)


    def __insertMesh(self, name, mesh_id, mesh_type, selection):
        """
        Add a new 'item' into the Hlist.
        """
        self.modelMesh.newData(name, mesh_id, mesh_type, selection)


    @pyqtSignature("")
    def slotAddMesh(self):
        """
        Add one monitoring point with these coordinates in the list in the Hlist
        The number of the monitoring point is added at the precedent one
        """
        mesh_id = self.mdl.addMesh()
        self.__insertMesh(self.mdl.getMeshLabel(mesh_id),
                                     mesh_id,
                                     self.mdl.getMeshType(mesh_id),
                                     self.mdl.getMeshLocation(mesh_id))


    @pyqtSignature("")
    def slotDeleteMesh(self):
        """
        Just delete the current selected entries from the Hlist and
        of course from the XML file.
        """
        list = []
        selectionModel = self.tableViewMesh.selectionModel()
        for index in selectionModel.selectedRows():
            mesh_id = self.modelMesh.getItem(index.row())['id']
            if int(mesh_id) < 0:
                title = self.tr("Warning")
                msg   = self.tr("You can't delete a default mesh")
                QMessageBox.information(self, title, msg)
                return
            list.append(str(mesh_id))

        self.mdl.deleteMesh(list)

        self.modelMesh.deleteAllData()
        list_mesh = []
        for mesh in self.mdl.getMeshIdList():
            if int(mesh)>0:
                list_mesh.append(mesh)
        new_id = 0
        for mesh in self.mdl.getMeshIdList():
            if int(mesh)<0:
                label = self.mdl.getMeshLabel(mesh)
                mesh_type = self.mdl.getMeshType(mesh)
                location = self.mdl.getMeshLocation(mesh)
                self.__insertMesh(label, mesh, mesh_type, location)
        for mesh in list_mesh:
            new_id = new_id + 1
            label = self.mdl.getMeshLabel(mesh)
            mesh_type = self.mdl.getMeshType(mesh)
            location = self.mdl.getMeshLocation(mesh)
            self.__insertMesh(label, str(new_id), mesh_type, location)


    @pyqtSignature("const QModelIndex&")
    def slotSelectMesh(self, index):
        #model = HeadLossesModel(self.case)
        cindex = self.tableViewMesh.currentIndex()
        if cindex != (-1,-1):
            row = cindex.row()
            mesh_id = self.modelMesh.getItem(row)['id']
            self.groupBoxVariable.show()
            self.groupBoxAssociatedWriter.show()
            if int(mesh_id) <0:
                self.checkBoxAllVariables.setEnabled(False)
                self.checkBoxAllVariables.setChecked(True)
                self.mdl.setMeshAllVariablesStatus(mesh_id,"on")
            else:
                self.checkBoxAllVariables.setEnabled(True)
                all_variables = self.mdl.getMeshAllVariablesStatus(mesh_id)
                if all_variables == 'on':
                    self.checkBoxAllVariables.setChecked(True)
                else :
                    self.checkBoxAllVariables.setChecked(False)
            self.showAssociatedWriterTable()


    @pyqtSignature("")
    def slotAllVariables(self):
        """
        Input INPDT0.
        """
        cindex = self.tableViewMesh.currentIndex()
        if cindex != (-1,-1):
            row = cindex.row()
            mesh_id = self.modelMesh.getItem(row)['id']
            if self.checkBoxAllVariables.isChecked():
                self.mdl.setMeshAllVariablesStatus(mesh_id, "on")
            else:
                self.mdl.setMeshAllVariablesStatus(mesh_id, "off")


    def __insertAssociatedWriter(self, name):
        """
        Add a new 'item' into the Hlist.
        """
        self.modelAssociatedWriter.newItem(name)


    @pyqtSignature("")
    def slotAddAssociatedWriter(self):
        """
        Add one monitoring point with these coordinates in the list in the Hlist
        The number of the monitoring point is added at the precedent one
        """
        cindex = self.tableViewMesh.currentIndex()
        if cindex != (-1,-1):
            row = cindex.row()
            mesh_id = self.modelMesh.getItem(row)['id']
            associated_writer_id = self.mdl.addAssociatedWriter(mesh_id)
            if associated_writer_id == None:
                title = self.tr("Warning")
                msg   = self.tr("Please create another writer\n"\
                                "before adding a new associated writer.")
                QMessageBox.information(self, title, msg)
                return
            self.__insertAssociatedWriter(self.mdl.getWriterLabel(associated_writer_id))


    @pyqtSignature("")
    def slotDeleteAssociatedWriter(self):
        """
        Just delete the current selected entries from the Hlist and
        of course from the XML file.
        """
        cindex = self.tableViewMesh.currentIndex()
        if cindex != (-1,-1):
            row = cindex.row()
            mesh_id = self.modelMesh.getItem(row)['id']
            selectionModel = self.tableViewAssociatedWriter.selectionModel()
            for index in selectionModel.selectedRows():
                writer_label = self.modelAssociatedWriter.getItem(index.row())
                writer_id = self.mdl.getWriterIdFromLabel(writer_label)
                self.mdl.deleteAssociatedWriter(mesh_id, writer_id)

                self.modelAssociatedWriter.deleteAllData()
                list_associated_writer = []
                for associated_writer in self.mdl.getAssociatedWriterIdList(mesh_id):
                    list_associated_writer.append(associated_writer)
                for associated_writer in list_associated_writer:
                    label = self.mdl.getWriterLabel(associated_writer)
                    self.__insertAssociatedWriter(label)


    @pyqtSignature("const QString &")
    def slotOutputProbeFmt(self, text):
        """
        INPUT choice of the output for the probes (.dat, .csv)
        """
        fmt = self.modelProbeFmt.dicoV2M[str(text)]
        log.debug("slotOutputProbeFmt-> fmt = %s" % fmt)
        self.mdl.setMonitoringPointFormat(fmt)


    @pyqtSignature("const QString &")
    def slotMonitoringPoint(self, text):
        """
        Input choice of the output of monitoring points files
        """
        histo = self.modelHisto.dicoV2M[str(text)]
        log.debug("slotMonitoringPoint-> histo = %s" % histo)
        self.mdl.setMonitoringPointType(histo)

        if histo == "None":
            nthist = -1
            self.mdl.setMonitoringPointFrequency(nthist)
            self.lineEditHisto.hide()
            self.lineEditFRHisto.hide()
            self.comboBoxProbeFmt.hide()
        else:
            self.comboBoxProbeFmt.show()

        if histo == "At each step":
            nthist = 1
            self.mdl.setMonitoringPointFrequency(nthist)
            self.lineEditHisto.hide()
            self.lineEditFRHisto.hide()

        if histo == "Frequency_h":
            self.lineEditHisto.show()
            self.lineEditHisto.setEnabled(True)
            nthist = self.mdl.getMonitoringPointFrequency()
            if nthist < 1:
                nthist = 1
                self.mdl.setMonitoringPointFrequency(nthist)
            self.lineEditHisto.setText(QString(str(nthist)))
            self.lineEditFRHisto.hide()

        if histo == "Frequency_h_x":
            self.lineEditHisto.hide()
            self.lineEditFRHisto.show()
            frlist = self.mdl.getMonitoringPointFrequencyTime()
            self.lineEditFRHisto.setText(QString(str(frlist)))


    @pyqtSignature("const QString &")
    def slotMonitoringPointFrequencyTime(self, text):
        """
        Input the frequency of the monitoring point output
        """
        n, ok = text.toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
            log.debug("slotMonitoringPointFrequencyTime-> FRHIST = %s" % n)
            self.mdl.setMonitoringPointFrequencyTime(n)


    @pyqtSignature("const QString &")
    def slotMonitoringPointFrequency(self, text):
        """
        Input the frequency of the monitoring point output
        """
        n, ok = text.toInt()
        if self.sender().validator().state == QValidator.Acceptable:
            log.debug("slotMonitoringPointFrequency-> NTHIST = %s" % n)
            self.mdl.setMonitoringPointFrequency(n)


    def __insertMonitoringPoint(self, num, X, Y, Z):
        """
        Add a new 'item' into the Hlist.
        """
        self.modelMonitoring.insertData(num, X, Y, Z)


    @pyqtSignature("")
    def slotAddMonitoringPoint(self):
        """
        Add one monitoring point with these coordinates in the list in the Hlist
        The number of the monitoring point is added at the precedent one
        """
        self.mdl.addMonitoringPoint(x=0.0, y=0.0, z=0.0)
        self.__insertMonitoringPoint(self.mdl.getNumberOfMonitoringPoints(),
                                     QString('0'),
                                     QString('0'),
                                     QString('0'))


    @pyqtSignature("")
    def slotDeleteMonitoringPoints(self):
        """
        Just delete the current selected entries from the Hlist and
        of course from the XML file.
        """
        list = []
        selectionModel = self.tableViewPoints.selectionModel()
        for index in selectionModel.selectedRows():
            name = index.row() + 1
            list.append(name)

        self.mdl.deleteMonitoringPoints(list)

        self.modelMonitoring.deleteAllData()
        for n in range(self.mdl.getNumberOfMonitoringPoints()):
            name = str(n+1)
            X, Y, Z = self.mdl.getMonitoringPointCoordinates(name)
            self.__insertMonitoringPoint(name, X, Y, Z)


    def isSteady(self):
        """
        """
        steady = 1
        from Pages.SteadyManagementModel import SteadyManagementModel

        if SteadyManagementModel(self.case).getSteadyFlowManagement() == 'on':
            steady = 0
        else:
            from Pages.TimeStepModel import TimeStepModel
            if TimeStepModel(self.case).getTimePassing() == 2:
                steady = 0

        return steady


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
