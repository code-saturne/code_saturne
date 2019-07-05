# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2019 EDF S.A.
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

import logging

#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------

from code_saturne.Base.QtCore    import *
from code_saturne.Base.QtGui     import *
from code_saturne.Base.QtWidgets import *
import os, re

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.model.Common import LABEL_LENGTH_MAX, GuiParam
from code_saturne.Base.QtPage import ComboModel, DoubleValidator, IntValidator
from code_saturne.Base.QtPage import RegExpValidator, to_qvariant
from code_saturne.Base.QtPage import from_qvariant, to_text_string
from code_saturne.Pages.OutputControlForm import Ui_OutputControlForm
from code_saturne.model.OutputControlModel import OutputControlModel
from code_saturne.Pages.QMeiEditorView import QMeiEditorView
from code_saturne.model.LagrangianModel import LagrangianModel
from code_saturne.model.NotebookModel import NotebookModel
from code_saturne.model.LocalizationModel import LocalizationModel

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
        rx = "[_A-Za-z0-9\(\)]{1," + str(LABEL_LENGTH_MAX-1) + "}"
        self.regExp = QRegExp(rx)
        v = RegExpValidator(editor, self.regExp)
        editor.setValidator(v)
        return editor


    def setEditorData(self, editor, index):
        editor.setAutoFillBackground(True)
        value = from_qvariant(index.model().data(index, Qt.DisplayRole), to_text_string)
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

                from code_saturne.Pages.VerifyExistenceLabelDialogView import VerifyExistenceLabelDialogView
                dialog = VerifyExistenceLabelDialogView(self.parent, default)
                if dialog.exec_():
                    result = dialog.get_result()
                    new_plabel = result['label']
                    log.debug("setModelData -> result = %s" % result)
                else:
                    new_plabel = self.old_plabel

            model.setData(index, to_qvariant(new_plabel), Qt.DisplayRole)


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
        rx = "[_A-Za-z0-9 \(\)]{1," + str(LABEL_LENGTH_MAX-1) + "}"
        self.regExp = QRegExp(rx)
        v = RegExpValidator(editor, self.regExp)
        editor.setValidator(v)
        return editor


    def setEditorData(self, editor, index):
        editor.setAutoFillBackground(True)
        value = from_qvariant(index.model().data(index, Qt.DisplayRole), to_text_string)
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

                from code_saturne.Pages.VerifyExistenceLabelDialogView import VerifyExistenceLabelDialogView
                dialog = VerifyExistenceLabelDialogView(self.parent, default)
                if dialog.exec_():
                    result = dialog.get_result()
                    new_plabel = result['label']
                    log.debug("setModelData -> result = %s" % result)
                else:
                    new_plabel = self.old_plabel

            model.setData(index, to_qvariant(new_plabel), Qt.DisplayRole)


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
        self.mdl = xml_model # TODO change this


    def createEditor(self, parent, option, index):
        editor = QComboBox(parent)
        editor.addItem("EnSight")
        editor.addItem("MED")
        editor.addItem("CGNS")
        editor.addItem("Catalyst")
        editor.addItem("Histogram")
        editor.addItem("Melissa")
        editor.addItem("CCM-IO")
        editor.installEventFilter(self)

        import cs_config
        cfg = cs_config.config()
        if cfg.libs['med'].have == "no":
            editor.setItemData(1, QColor(Qt.red), Qt.TextColorRole);
        if cfg.libs['cgns'].have == "no":
            editor.setItemData(2, QColor(Qt.red), Qt.TextColorRole);
        if cfg.libs['catalyst'].have == "no":
            editor.setItemData(3, QColor(Qt.red), Qt.TextColorRole);
        if cfg.libs['melissa'].have == "no":
            editor.setItemData(5, QColor(Qt.red), Qt.TextColorRole);
        if cfg.libs['ccm'].have == "no":
            editor.setItemData(6, QColor(Qt.red), Qt.TextColorRole);
        return editor


    def setEditorData(self, comboBox, index):
        dico = {"ensight": 0, "med": 1, "cgns": 2, "catalyst": 3,
                "histogram": 4, "melissa": 5, "ccm": 6}
        row = index.row()
        string = index.model().dataWriter[row]['format']
        idx = dico[string]
        comboBox.setCurrentIndex(idx)


    def setModelData(self, comboBox, model, index):
        value = comboBox.currentText()
        selectionModel = self.parent.selectionModel()
        for idx in selectionModel.selectedIndexes():
            if idx.column() == index.column():
                model.setData(idx, to_qvariant(value))


#-------------------------------------------------------------------------------
# Combo box delegate for the writer format
#-------------------------------------------------------------------------------

class TypeMeshDelegate(QItemDelegate):
    """
    Use of a combo box in the table.
    """
    def __init__(self, parent=None, xml_model=None, typ=0):
        super(TypeMeshDelegate, self).__init__(parent)
        self.parent = parent
        self.mdl = xml_model # TODO review this
        self.lag = typ


    def createEditor(self, parent, option, index):
        editor = QComboBox(parent)
        if self.lag == 0:
            editor.addItem("cells")
            editor.addItem("interior faces")
            editor.addItem("boundary faces")
            editor.addItem("volume zone")
            editor.addItem("boundary zone")

            vlm = LocalizationModel('VolumicZone', self.mdl.case)
            if len(vlm.getZones()) == 0:
                idx = editor.findText('volume zone')
                editor.model().item(idx).setEnabled(False)
            blm = LocalizationModel('BoundaryZone', self.mdl.case)
            if len(blm.getZones()) == 0:
                idx = editor.findText('boundary zone')
                editor.model().item(idx).setEnabled(False)

        else:
            editor.addItem("particles")
            editor.addItem("trajectories")
        editor.installEventFilter(self)
        return editor


    def setEditorData(self, comboBox, index):
        if self.lag == 0:
            dico = {"cells": 0, "interior_faces": 1, "boundary_faces": 2,
                    "VolumicZone": 3, "BoundaryZone": 4}
        else:
            dico = {"particles": 0, "trajectories": 1}
        row = index.row()
        string = index.model().dataMesh[row]['type']
        idx = dico[string]
        comboBox.setCurrentIndex(idx)


    def setModelData(self, comboBox, model, index):
        value = comboBox.currentText()
        selectionModel = self.parent.selectionModel()
        for idx in selectionModel.selectedIndexes():
            if idx.column() == index.column():
                model.setData(idx, to_qvariant(value))


    def paint(self, painter, option, index):
        row = index.row()
        meshtype = index.model().dataMesh[row]['type']
        isValid = meshtype != None and meshtype != ''

        if isValid:
            QItemDelegate.paint(self, painter, option, index)

        else:
            painter.save()
            # set background color
            if option.state & QStyle.State_Selected:
                painter.setBrush(QBrush(Qt.darkRed))
            else:
                painter.setBrush(QBrush(Qt.red))
            # set text color
            painter.setPen(QPen(Qt.NoPen))
            painter.drawRect(option.rect)
            painter.setPen(QPen(Qt.black))
            value = index.data(Qt.DisplayRole)
            if value.isValid():
                text = from_qvariant(value, to_text_string)
                painter.drawText(option.rect, Qt.AlignLeft, text)
            painter.restore()


#-------------------------------------------------------------------------------
# QLineEdit delegate for location
#-------------------------------------------------------------------------------

class LocationSelectorDelegate(QItemDelegate):
    def __init__(self, parent, mdl):
        super(LocationSelectorDelegate, self).__init__(parent)
        self.parent = parent
        self.mdl = mdl
        self.editor_type = 'line'


    def createEditor(self, parent, option, index):

        row = index.row()
        mesh_id = self.mdl.getMeshIdList()[row]
        mesh_type = self.mdl.getMeshType(mesh_id)

        if mesh_type in ['BoundaryZone', 'VolumicZone']:
            editor = QComboBox(parent)
            lm = LocalizationModel(mesh_type, self.mdl.case)
            for zone in lm.getZones():
                editor.addItem(zone.getLabel())

            self.editor_type = 'combo'
        else:
            editor = QLineEdit(parent)

            # Autocompletion for selection criteria!
            comp_list = ['all[]']
            comp_list += ['plane[a, b, c, d, epsilon]',
                          'plane[a, b, c, d, inside]',
                          'plane[a, b, c, d, outside]',
                          'plane[n_x, n_y, n_z, x0, y0, z0, epsilon]',
                          'plane[n_x, n_y, n_z, x0, y0, z0, inside]',
                          'plane[n_x, n_y, n_z, x0, y0, z0, outside]',
                          'box[xmin, ymin, zmin, xmax, ymax, zmax]',
                          'box[x0, y0, z0, dx1, dy1, dz1, dx2, dy2, dz2, dx3, dy3, dz3]',
                          'cylinder[x0, y0, z0, x1, y1, z1, radius]',
                          'sphere[x_c, y_c, z_c, radius]']

            completer = QCompleter()
            editor.setCompleter(completer)
            model = QStringListModel()
            completer.setModel(model)
            model.setStringList(comp_list)

            self.editor_type = 'line'

        return editor


    def setEditorData(self, editor, index):
        # This line is used to avoid an overlay of old and new text
        editor.setAutoFillBackground(True)
        self.value = from_qvariant(index.model().data(index, Qt.DisplayRole),
                                   to_text_string)
        if self.editor_type == 'line':
            editor.setText(self.value)
        else:
            idx = editor.findText(self.value)
            editor.setCurrentIndex(idx)


    def setModelData(self, editor, model, index):

        if self.editor_type == 'line':
            value = editor.text()
        else:
            value = editor.currentText()

        if str(value) == "" :
           title = self.tr("Warning")
           msg   = self.tr("Please give a location")
           QMessageBox.information(self.parent, title, msg)
           return

        if str(value) != "" :
            model.setData(index, to_qvariant(value), Qt.DisplayRole)


#-------------------------------------------------------------------------------
# QLineEdit delegate for density
#-------------------------------------------------------------------------------

class DensitySelectorDelegate(QItemDelegate):
    def __init__(self, parent, mdl):
        super(DensitySelectorDelegate, self).__init__(parent)
        self.parent = parent
        self.mdl = mdl


    def createEditor(self, parent, option, index):
        editor = QLineEdit(parent)
        validator = DoubleValidator(editor, min=0., max=1.)
        editor.setValidator(validator)
        return editor


    def setEditorData(self, editor, index):
        editor.setAutoFillBackground(True)
        self.value = from_qvariant(index.model().data(index, Qt.DisplayRole), to_text_string)
        editor.setText(self.value)


    def setModelData(self, editor, model, index):
        value = editor.text()

        if str(value) == "" :
           title = self.tr("Warning")
           msg   = self.tr("Please give a density")
           QMessageBox.information(self.parent, title, msg)
           return

        if str(value) != "" :
            model.setData(index, to_qvariant(value), Qt.DisplayRole)


#-------------------------------------------------------------------------------
# Combo box delegate for the variance
#-------------------------------------------------------------------------------

class AssociatedWriterDelegate(QItemDelegate):
    """
    Use of a combo box in the table.
    """
    def __init__(self, parent, lag = 0):
        super(AssociatedWriterDelegate, self).__init__(parent)
        self.parent   = parent
        self.lagrangian = lag


    def createEditor(self, parent, option, index):
        editor = QComboBox(parent)
        self.modelCombo = ComboModel(editor, 1, 1)
        editor.installEventFilter(self)
        return editor


    def setEditorData(self, editor, index):
        editor.setAutoFillBackground(True)
        l1 = index.model().mdl.getWriterLabelList()
        for s in l1:
            idx = index.model().mdl.getWriterIdFromLabel(s)
            if self.lagrangian == 0:
                if int(idx) > -2:
                    self.modelCombo.addItem(s, s)
            else:
                if int(idx) > 0 or int(idx) == -3 or int(idx) == -4:
                    self.modelCombo.addItem(s, s)


    def setModelData(self, comboBox, model, index):
        txt = str(comboBox.currentText())
        value = self.modelCombo.dicoV2M[txt]
        model.setData(index, to_qvariant(value), Qt.DisplayRole)


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
        self.dicoV2M= {"cells": 'cells',
                       "interior faces" : 'interior_faces',
                       "boundary faces": 'boundary_faces',
                       "volume zone" : 'VolumicZone',
                       "boundary zone" : 'BoundaryZone'}
        self.dicoM2V= {"cells" : 'cells',
                       "interior_faces" : 'interior faces',
                       "boundary_faces": 'boundary faces',
                       "VolumicZone": 'volume zone',
                       "BoundaryZone": 'boundary zone'}
        type_list = ["cells", "interior_faces", "boundary_faces",
                     "VolumicZone", "BoundaryZone"]
        for id in self.mdl.getMeshIdList():
            dico  = {}
            dico['name'] = self.mdl.getMeshLabel(id)
            dico['id'] = id
            dico['type'] = self.mdl.getMeshType(id)
            dico['location'] = self.mdl.getMeshLocation(id)
            if dico['type'] in type_list:
                row = self.rowCount()
                self.setRowCount(row + 1)
                self.dataMesh.append(dico)
                if int(id) < 0:
                    self.defaultItem.append(row)
            log.debug("populateModel-> dataSolver = %s" % dico)


    def data(self, index, role):
        if not index.isValid():
            return to_qvariant()

        if role == Qt.DisplayRole:

            row = index.row()
            col = index.column()
            dico = self.dataMesh[row]

            if index.column() == 0:
                return to_qvariant(dico['name'])
            elif index.column() == 1:
                return to_qvariant(dico['id'])
            elif index.column() == 2:
                return to_qvariant(self.dicoM2V[dico['type']])
            elif index.column() == 3:
                return to_qvariant(dico['location'])
            else:
                return to_qvariant()

        elif role == Qt.TextAlignmentRole:
            if index.column() != 3:
                return to_qvariant(Qt.AlignCenter)
            else:
                return to_qvariant(Qt.AlignLeft | Qt.AlignVCenter)

        return to_qvariant()


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled
        # default item
        col_id = index.column()
        if index.row() in self.defaultItem:
            if col_id == 1 or col_id == 2:
                return Qt.ItemIsEnabled | Qt.ItemIsSelectable
        else:
            if col_id == 1:
                return Qt.ItemIsEnabled | Qt.ItemIsSelectable
        return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable


    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            if section == 0:
                return to_qvariant(self.tr("Name"))
            elif section == 1:
                return to_qvariant(self.tr("Id"))
            elif section == 2:
                return to_qvariant(self.tr("Type"))
            elif section == 3:
                return to_qvariant(self.tr("Selection Criteria"))
        return to_qvariant()


    def setData(self, index, value, role=None):

        # Update the row in the table
        row = index.row()
        col = index.column()

        # Lists to better handle the difference between elements and zones
        elts_list = ['cells', 'interior_faces', 'boundary_faces']
        zone_list = ['VolumicZone', 'BoundaryZone']

        # Label
        if col == 0:
            old_plabel = self.dataMesh[row]['name']
            new_plabel = str(from_qvariant(value, to_text_string))
            self.dataMesh[row]['name'] = new_plabel
            self.mdl.setMeshLabel(str(self.dataMesh[row]['id']), new_plabel)

        if index.column() == 2:
            old_type = self.dataMesh[row]['type']
            new_type = self.dicoV2M[str(from_qvariant(value, to_text_string))]

            if new_type in zone_list and old_type != new_type:
                self.dataMesh[row]['location'] = \
                    LocalizationModel(new_type, self.mdl.case).getZones()[0].getLabel()

            elif new_type in elts_list and old_type in zone_list:
                self.dataMesh[row]['location'] = 'all[]'

            self.dataMesh[row]['type'] = new_type
            self.mdl.setMeshType(self.dataMesh[row]['id'], self.dataMesh[row]['type'])

        if index.column() == 3:
            new_location = str(from_qvariant(value, to_text_string))
            self.dataMesh[row]['location'] = new_location
            self.mdl.setMeshLocation(self.dataMesh[row]['id'], new_location)
        self.dataChanged.emit(index, index)
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
# QStandardItemModel for Lagrangian Mesh QTableView
#-------------------------------------------------------------------------------

class StandardItemModelLagrangianMesh(QStandardItemModel):
    def __init__(self, mdl):
        """
        """
        QStandardItemModel.__init__(self)

        self.setColumnCount(5)
        self.dataMesh = []
        self.mdl = mdl
        self.defaultItem = []
        self.populateModel()


    def populateModel(self):

        self.dicoV2M= {"particles": 'particles',
                       "trajectories": 'trajectories'}
        self.dicoM2V= {"particles" : 'particles',
                       "trajectories": 'trajectories'}
        for idx in self.mdl.getMeshIdList():
            dico  = {}
            dico['name'] = self.mdl.getMeshLabel(idx)
            dico['id'] = idx
            dico['type'] = self.mdl.getLagrangianMeshType(idx)
            dico['density'] = self.mdl.getMeshDensity(idx)
            dico['location'] = self.mdl.getMeshLocation(idx)

            if dico['type'] in ['particles', 'trajectories']:
                row = self.rowCount()
                self.setRowCount(row + 1)
                self.dataMesh.append(dico)
                if int(idx) < 0:
                    self.defaultItem.append(row)
            log.debug("populateModel-> dataSolver = %s" % dico)


    def data(self, index, role):
        if not index.isValid():
            return to_qvariant()

        if role == Qt.DisplayRole:

            row = index.row()
            col = index.column()
            dico = self.dataMesh[row]

            if index.column() == 0:
                return to_qvariant(dico['name'])
            elif index.column() == 1:
                return to_qvariant(dico['id'])
            elif index.column() == 2:
                return to_qvariant(self.dicoM2V[dico['type']])
            elif index.column() == 3:
                return to_qvariant(dico['density'])
            elif index.column() == 4:
                return to_qvariant(dico['location'])
            else:
                return to_qvariant()

        elif role == Qt.TextAlignmentRole:
            if index.column() != 4:
                return to_qvariant(Qt.AlignCenter)
            else:
                return to_qvariant(Qt.AlignLeft | Qt.AlignVCenter)

        return to_qvariant()


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled
        # default item
        col_id = index.column()
        if index.row() in self.defaultItem:
            if col_id == 1 or col_id == 2:
                return Qt.ItemIsEnabled | Qt.ItemIsSelectable
        else:
            if col_id == 1:
                return Qt.ItemIsEnabled | Qt.ItemIsSelectable
        return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable


    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            if section == 0:
                return to_qvariant(self.tr("Name"))
            elif section == 1:
                return to_qvariant(self.tr("Id"))
            elif section == 2:
                return to_qvariant(self.tr("Type"))
            elif section == 3:
                return to_qvariant(self.tr("Density"))
            elif section == 4:
                return to_qvariant(self.tr("Selection Criteria"))
        return to_qvariant()


    def setData(self, index, value, role=None):

        # Update the row in the table
        row = index.row()
        col = index.column()

        # Label
        if col == 0:
            old_plabel = self.dataMesh[row]['name']
            new_plabel = str(from_qvariant(value, to_text_string))
            self.dataMesh[row]['name'] = new_plabel
            self.mdl.setMeshLabel(str(self.dataMesh[row]['id']), new_plabel)

        if index.column() == 2:
            self.dataMesh[row]['type'] = self.dicoV2M[str(from_qvariant(value, to_text_string))]
            self.mdl.setLagrangianMeshType(self.dataMesh[row]['id'], self.dataMesh[row]['type'])

        if index.column() == 3:
            self.dataMesh[row]['density'] = str(from_qvariant(value, to_text_string))
            self.mdl.setMeshDensity(self.dataMesh[row]['id'], self.dataMesh[row]['density'])

        if index.column() == 4:
            new_location = str(from_qvariant(value, to_text_string))
            self.dataMesh[row]['location'] = new_location
            self.mdl.setMeshLocation(self.dataMesh[row]['id'], new_location)

        self.dataChanged.emit(index, index)
        return True


    def newData(self, name, mesh_id, mesh_type, density, location):
        """
        Add a new 'item' into the table.
        """
        dico = {}
        dico['name'] = name
        dico['id'] = mesh_id
        dico['type'] = mesh_type
        dico['density'] = density
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
    def __init__(self, mdl, parent):
        """
        """
        QStandardItemModel.__init__(self)

        self.setColumnCount(4)
        self.dataWriter = []
        self.mdl = mdl
        self.defaultItem = []
        self.parent = parent
        self.populateModel()

    def populateModel(self):
        self.dicoV2M= {"EnSight": 'ensight',
                       "MED" : 'med',
                       "CGNS": 'cgns',
                       "Catalyst": 'catalyst',
                       "Histogram": 'histogram',
                       "Melissa": 'melissa',
                       "CCM-IO": 'ccm'}
        self.dicoM2V= {"ensight" : 'EnSight',
                       "med" : 'MED',
                       "cgns": 'CGNS',
                       "catalyst": 'Catalyst',
                       "histogram": 'Histogram',
                       "melissa": 'Melissa',
                       "ccm": 'CCM-IO'}
        for id in self.mdl.getWriterIdList():
            row = self.rowCount()
            self.setRowCount(row + 1)

            dico = {}
            dico['name'] = self.mdl.getWriterLabel(id)
            dico['id'] = id
            dico['format'] = self.mdl.getWriterFormat(id)
            dico['directory'] = self.mdl.getWriterDirectory(id)

            self.dataWriter.append(dico)
            if int(id)<0:
                self.defaultItem.append(row)
            log.debug("populateModel-> dataSolver = %s" % dico)

    def data(self, index, role):
        if not index.isValid():
            return to_qvariant()

        if role == Qt.DisplayRole:

            row = index.row()
            col = index.column()
            dico = self.dataWriter[row]

            if index.column() == 0:
                return to_qvariant(dico['name'])
            elif index.column() == 1:
                return to_qvariant(dico['id'])
            elif index.column() == 2:
                return to_qvariant(self.dicoM2V[dico['format']])
            elif index.column() == 3:
                return to_qvariant(dico['directory'])
            else:
                return to_qvariant()

        elif role == Qt.TextAlignmentRole:
            return to_qvariant(Qt.AlignCenter)

        return to_qvariant()


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled
        # default item
        if index.column() == 1:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable
        else:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable


    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            if section == 0:
                return to_qvariant(self.tr("Name"))
            elif section == 1:
                return to_qvariant(self.tr("Id"))
            elif section == 2:
                return to_qvariant(self.tr("Format"))
            elif section == 3:
                return to_qvariant(self.tr("Directory"))
        return to_qvariant()


    def setData(self, index, value, role=None):

        # Update the row in the table
        row = index.row()
        col = index.column()

        writer_id = self.dataWriter[row]['id']
        # Label
        if col == 0:
            old_plabel = self.dataWriter[row]['name']
            new_plabel = str(from_qvariant(value, to_text_string))
            self.dataWriter[row]['name'] = new_plabel
            self.mdl.setWriterLabel(writer_id, new_plabel)

        elif col == 2:
            f_old = self.mdl.getWriterFormat(writer_id)
            self.dataWriter[row]['format'] = self.dicoV2M[str(from_qvariant(value, to_text_string))]
            if self.dataWriter[row]['format'] != f_old:
                self.mdl.setWriterFormat(writer_id,
                                         self.dataWriter[row]['format'])
                options = self.mdl.getWriterOptions(writer_id)
                new_options = ""
                if 'separate_meshes' in options:
                  new_options = "separate_meshes"
                self.mdl.setWriterOptions(writer_id, new_options)
        elif col == 3:
            old_rep = self.dataWriter[row]['directory']
            new_rep = str(from_qvariant(value, to_text_string))
            self.dataWriter[row]['directory'] = new_rep
            self.mdl.setWriterDirectory(writer_id, new_rep)


        self.dataChanged.emit(index, index)
        return True


    def newData(self, name, writer_id, writer_format, directory):
        """
        Add a new 'item' into the table.
        """
        dico = {}
        dico['name'] = name
        dico['id'] = writer_id
        dico['format'] = writer_format
        dico['directory'] = directory
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
            return to_qvariant()

        row = index.row()
        col = index.column()

        if role == Qt.DisplayRole:
            return to_qvariant(self._data[row])

        return to_qvariant()


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled
        return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable


    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return to_qvariant(self.headers[section])
        return to_qvariant()


    def setData(self, index, value, role):
        if not index.isValid():
            return Qt.ItemIsEnabled
        # Update the row in the table
        row = index.row()
        col = index.column()

        # Label
        if col == 0:
            writer = str(from_qvariant(value, to_text_string))
            self._data[row] = writer
            writer_list = []
            for r in range(self.rowCount()):
                writer_list.append(self.mdl.getWriterIdFromLabel(self._data[r]))
            self.mdl.setAssociatedWriterChoice(self.mesh_id, writer_list)

        self.dataChanged.emit(index, index)
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
            return to_qvariant()

        if role == Qt.DisplayRole:

            row = index.row()
            dico = self.dataMonitoring[row]

            if index.column() == 0:
                return to_qvariant(dico['n'])
            elif index.column() == 1:
                return to_qvariant(dico['X'])
            elif index.column() == 2:
                return to_qvariant(dico['Y'])
            elif index.column() == 3:
                return to_qvariant(dico['Z'])
            else:
                return to_qvariant()

        elif role == Qt.TextAlignmentRole:
            return to_qvariant(Qt.AlignCenter)

        return to_qvariant()


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
                return to_qvariant(self.tr("n"))
            elif section == 1:
                return to_qvariant(self.tr("X"))
            elif section == 2:
                return to_qvariant(self.tr("Y"))
            elif section == 3:
                return to_qvariant(self.tr("Z"))
        return to_qvariant()


    def setData(self, index, value, role=None):
        row = index.row()
        if index.column() == 0:
            n = from_qvariant(value, int)
            self.dataMonitoring[row]['n'] = n
        elif index.column() == 1:
            X = from_qvariant(value, float)
            self.dataMonitoring[row]['X'] = X
        elif index.column() == 2:
            Y = from_qvariant(value, float)
            self.dataMonitoring[row]['Y'] = Y
        elif index.column() == 3:
            Z = from_qvariant(value, float)
            self.dataMonitoring[row]['Z'] = Z

        self.dataChanged.emit(index, index)
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


    def replaceData(self, row, label, X, Y, Z):
        """
        Replace value in an existing 'item' into the table.
        """
        self.dataMonitoring[row]['n'] = label
        self.dataMonitoring[row]['X'] = X
        self.dataMonitoring[row]['Y'] = Y
        self.dataMonitoring[row]['Z'] = Z

        row = self.rowCount()
        self.setRowCount(row)


    def getLabel(self, index):
        return self.getData(index)['n']


    def getData(self, index):
        row = index.row()
        dico = self.dataMonitoring[row]
        return dico


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
    def __init__(self, parent=None, case=None, model=None):
        """ Construtor.

        @param: parent ancestor object
        @model: monitoring points model
        """
        super(MonitoringPointDelegate, self).__init__(parent)
        self.table = parent
        self.case = case
        self.mdl = model


    def createEditor(self, parent, option, index):
        if index.column() == 0:
            editor = QFrame(parent)
        else:
            editor = QLineEdit(parent)
            editor.setValidator(DoubleValidator(editor))
            editor.setFrame(False)
            editor.returnPressed.connect(self.commitAndCloseEditor)
            editor.setCursorPosition(0)
        return editor


    def commitAndCloseEditor(self):
        editor = self.sender()
        if isinstance(editor, QLineEdit):
            self.commitData.emit(editor)
            self.closeEditor.emit(editor, QAbstractItemDelegate.NoHint)


    def setEditorData(self, editor, index):
        editor.setAutoFillBackground(True)
        text = from_qvariant(index.model().data(index, Qt.DisplayRole), to_text_string)
        if isinstance(editor, QLineEdit):
            editor.setText(text)


    def setModelData(self, editor, model, index):
        if isinstance(editor, QLineEdit):
            if not editor.isModified():
                return

            item = editor.text()
            selectionModel = self.table.selectionModel()
            for index in selectionModel.selectedRows(index.column()):
                model.setData(index, to_qvariant(item), Qt.DisplayRole)
                dico = model.dataMonitoring[index.row()]
                x = float(dico['X'])
                y = float(dico['Y'])
                z = float(dico['Z'])
                label = str(dico['n'])
                self.mdl.replaceMonitoringPointCoordinates(label, x, y, z)
                if self.case['probes']:
                    self.case['probes'].updateLocation(label, [x, y, z])

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
        self.case.undoStopGlobal()
        self.mdl = OutputControlModel(self.case)
        self.notebook = NotebookModel(self.case)

        import cs_config
        cfg = cs_config.config()
        no_catalyst = False
        if cfg.libs['catalyst'].have == "no":
            no_catalyst = True

        if self.case['run_type'] == 'standard':
            # lagrangian model
            self.lag_mdl = LagrangianModel(self.case)

        # Combo models

        self.modelOutput         = ComboModel(self.comboBoxOutput,3,1)
        self.modelNTLAL          = ComboModel(self.comboBoxNTLAL,3,1)
        self.modelFrequency      = ComboModel(self.comboBoxFrequency,4,1)
        self.modelTimeDependency = ComboModel(self.comboBoxTimeDependency,3,1)
        self.modelFormatE        = ComboModel(self.comboBoxFormatE,2,1)
        self.modelFormatH        = ComboModel(self.comboBoxFormatH,2,1)
        self.modelPolygon        = ComboModel(self.comboBoxPolygon,3,1)
        self.modelPolyhedra      = ComboModel(self.comboBoxPolyhedra,3,1)
        self.modelTimePlot       = ComboModel(self.comboBoxTimePlot,3,1)
        self.modelProbeFmt       = ComboModel(self.comboBoxProbeFmt,2,1)

        self.modelOutput.addItem(self.tr("No output"), 'None')
        self.modelOutput.addItem(self.tr("Output listing at each time step"), 'At each step')
        self.modelOutput.addItem(self.tr("Output every 'n' time steps"), 'Frequency_l')

        self.modelNTLAL.addItem(self.tr("No output"), 'None')
        self.modelNTLAL.addItem(self.tr("Output listing at each time step"), 'At each step')
        self.modelNTLAL.addItem(self.tr("Output every 'n' time steps"), 'Frequency_l')

        self.modelFrequency.addItem(self.tr("No periodic output"), 'none')
        self.modelFrequency.addItem(self.tr("Output every 'n' time steps"), 'time_step')
        self.modelFrequency.addItem(self.tr("Output every 'x' seconds"), 'time_value')
        self.modelFrequency.addItem(self.tr("Output using a formula"), 'formula')

        self.modelTimeDependency.addItem(self.tr("Fixed mesh"), 'fixed_mesh')
        self.modelTimeDependency.addItem(self.tr("Transient coordinates"), 'transient_coordinates')
        self.modelTimeDependency.addItem(self.tr("Transient connectivity"), 'transient_connectivity')

        self.modelFormatE.addItem(self.tr("binary (native)"), 'binary')
        self.modelFormatE.addItem(self.tr("binary (big-endian)"), 'big_endian')
        self.modelFormatE.addItem(self.tr("text"), 'text')

        self.modelFormatH.addItem(self.tr("text"), 'text')
        self.modelFormatH.addItem(self.tr("TeX (TikZ)"), 'tex')
        self.modelFormatH.addItem(self.tr("PNG"), 'png', warn=no_catalyst)

        self.modelPolygon.addItem(self.tr("display"), 'display')
        self.modelPolygon.addItem(self.tr("discard"), 'discard_polygons')
        self.modelPolygon.addItem(self.tr("subdivide"), 'divide_polygons')

        self.modelPolyhedra.addItem(self.tr("display"), 'display')
        self.modelPolyhedra.addItem(self.tr("discard"), 'discard_polyhedra')
        self.modelPolyhedra.addItem(self.tr("subdivide"), 'divide_polyhedra')

        self.modelTimePlot.addItem(self.tr("No monitoring points file"), 'None')
        self.modelTimePlot.addItem(self.tr("Monitoring points files at each time step"), 'At each step')
        self.modelTimePlot.addItem(self.tr("Monitoring points file every 'n' time steps"), 'Frequency_h')
        self.modelTimePlot.addItem(self.tr("Monitoring points file every 'x' time_value(s)"), 'Frequency_h_x')

        self.modelProbeFmt.addItem(self.tr(".dat"), 'DAT')
        self.modelProbeFmt.addItem(self.tr(".csv"), 'CSV')

        # Hide time frequency (in s) when calculation is steady
        if self.case['run_type'] == 'standard':
            if self.isSteady() != 0:
                self.modelTimePlot.disableItem(3)

        # Model for QTableView

        self.modelMonitoring = StandardItemModelMonitoring()
        self.tableViewPoints.setModel(self.modelMonitoring)
        self.tableViewPoints.resizeColumnToContents(0)
        self.tableViewPoints.resizeRowsToContents()
        self.tableViewPoints.setAlternatingRowColors(True)
        self.tableViewPoints.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.tableViewPoints.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.tableViewPoints.setEditTriggers(QAbstractItemView.DoubleClicked)
        if QT_API == "PYQT4":
            self.tableViewPoints.horizontalHeader().setResizeMode(QHeaderView.Stretch)
        elif QT_API == "PYQT5":
            self.tableViewPoints.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.delegate = MonitoringPointDelegate(self.tableViewPoints, self.case, self.mdl)
        self.tableViewPoints.setItemDelegate(self.delegate)

        self.modelWriter = StandardItemModelWriter(self.mdl, parent)
        self.tableViewWriter.setModel(self.modelWriter)
        self.tableViewWriter.resizeColumnToContents(0)
        self.tableViewWriter.resizeRowsToContents()
        self.tableViewWriter.setAlternatingRowColors(True)
        self.tableViewWriter.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.tableViewWriter.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.tableViewWriter.setEditTriggers(QAbstractItemView.DoubleClicked)
        if QT_API == "PYQT4":
            self.tableViewWriter.horizontalHeader().setResizeMode(QHeaderView.Stretch)
        elif QT_API == "PYQT5":
            self.tableViewWriter.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)

        delegate_label_writer = LabelWriterDelegate(self.tableViewWriter)
        self.tableViewWriter.setItemDelegateForColumn(0, delegate_label_writer)
        self.tableViewWriter.setItemDelegateForColumn(3, delegate_label_writer)
        delegate_format = FormatWriterDelegate(self.tableViewWriter, self.mdl)
        self.tableViewWriter.setItemDelegateForColumn(2, delegate_format)

        # mesh tab
        self.modelMesh = StandardItemModelMesh(self.mdl)
        self.tableViewMesh.setModel(self.modelMesh)
        self.tableViewMesh.resizeColumnToContents(0)
        self.tableViewMesh.resizeRowsToContents()
        self.tableViewMesh.setAlternatingRowColors(True)
        self.tableViewMesh.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.tableViewMesh.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.tableViewMesh.setEditTriggers(QAbstractItemView.DoubleClicked)
        if QT_API == "PYQT4":
            self.tableViewMesh.horizontalHeader().setResizeMode(0, QHeaderView.ResizeToContents)
            self.tableViewMesh.horizontalHeader().setResizeMode(1, QHeaderView.ResizeToContents)
            self.tableViewMesh.horizontalHeader().setResizeMode(2, QHeaderView.ResizeToContents)
            self.tableViewMesh.horizontalHeader().setResizeMode(3, QHeaderView.ResizeToContents)
        elif QT_API == "PYQT5":
            self.tableViewMesh.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeToContents)
            self.tableViewMesh.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeToContents)
            self.tableViewMesh.horizontalHeader().setSectionResizeMode(2, QHeaderView.ResizeToContents)
            self.tableViewMesh.horizontalHeader().setSectionResizeMode(3, QHeaderView.ResizeToContents)
        self.tableViewMesh.horizontalHeader().setStretchLastSection(True)

        delegate_label_mesh = LabelMeshDelegate(self.tableViewMesh)
        self.tableViewMesh.setItemDelegateForColumn(0, delegate_label_mesh)
        delegate_type = TypeMeshDelegate(self.tableViewMesh, self.mdl, 0)
        self.tableViewMesh.setItemDelegateForColumn(2, delegate_type)
        delegate_location = LocationSelectorDelegate(self.tableViewMesh, self.mdl)
        self.tableViewMesh.setItemDelegateForColumn(3, delegate_location)

        # lagrangian mesh tab
        self.modelLagrangianMesh = StandardItemModelLagrangianMesh(self.mdl)
        self.tableViewLagrangianMesh.setModel(self.modelLagrangianMesh)
        self.tableViewLagrangianMesh.resizeColumnToContents(0)
        self.tableViewLagrangianMesh.resizeRowsToContents()
        self.tableViewLagrangianMesh.setAlternatingRowColors(True)
        self.tableViewLagrangianMesh.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.tableViewLagrangianMesh.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.tableViewLagrangianMesh.setEditTriggers(QAbstractItemView.DoubleClicked)
        if QT_API == "PYQT4":
            self.tableViewLagrangianMesh.horizontalHeader().setResizeMode(0, QHeaderView.ResizeToContents)
            self.tableViewLagrangianMesh.horizontalHeader().setResizeMode(1, QHeaderView.ResizeToContents)
            self.tableViewLagrangianMesh.horizontalHeader().setResizeMode(2, QHeaderView.ResizeToContents)
        elif QT_API == "PYQT5":
            self.tableViewLagrangianMesh.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeToContents)
            self.tableViewLagrangianMesh.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeToContents)
            self.tableViewLagrangianMesh.horizontalHeader().setSectionResizeMode(2, QHeaderView.ResizeToContents)
        self.tableViewLagrangianMesh.horizontalHeader().setStretchLastSection(True)

        delegate_label_lag_mesh = LabelMeshDelegate(self.tableViewLagrangianMesh)
        self.tableViewLagrangianMesh.setItemDelegateForColumn(0, delegate_label_lag_mesh)
        delegate_lag_type = TypeMeshDelegate(self.tableViewLagrangianMesh, self.mdl, 1)
        self.tableViewLagrangianMesh.setItemDelegateForColumn(2, delegate_lag_type)
        delegate_lag_density = DensitySelectorDelegate(self.tableViewLagrangianMesh, self.mdl)
        self.tableViewLagrangianMesh.setItemDelegateForColumn(3, delegate_lag_density)
        delegate_lag_location = LocationSelectorDelegate(self.tableViewLagrangianMesh, self.mdl)
        self.tableViewLagrangianMesh.setItemDelegateForColumn(4, delegate_lag_location)

        # Connections

        self.tabWidget.currentChanged[int].connect(self.slotchanged)
        self.modelWriter.dataChanged.connect(self.slotDataChanged)
        self.tableViewMesh.clicked.connect(self.slotSelectMesh)
        self.tableViewLagrangianMesh.clicked.connect(self.slotSelectLagrangianMesh)
        self.tableViewWriter.clicked.connect(self.slotSelectWriter)
        self.checkBoxOutputStart.clicked.connect(self.slotWriterOutputStart)
        self.checkBoxOutputEnd.clicked.connect(self.slotWriterOutputEnd)
        self.checkBoxSeparateMeshes.clicked.connect(self.slotWriterSeparateMeshes)
        self.checkBoxAllowParallelIO.clicked.connect(self.slotWriterParallelIO)
        self.checkBoxAllVariables.clicked.connect(self.slotAllVariables)
        self.checkBoxAllLagrangianVariables.clicked.connect(self.slotAllLagrangianVariables)
        self.pushButtonFrequency.clicked.connect(self.slotWriterFrequencyFormula)

        self.pushButtonAddWriter.clicked.connect(self.slotAddWriter)
        self.pushButtonDeleteWriter.clicked.connect(self.slotDeleteWriter)
        self.pushButtonAddMesh.clicked.connect(self.slotAddMesh)
        self.pushButtonDeleteMesh.clicked.connect(self.slotDeleteMesh)
        self.pushButtonAddAssociatedWriter.clicked.connect(self.slotAddAssociatedWriter)
        self.pushButtonDeleteAssociatedWriter.clicked.connect(self.slotDeleteAssociatedWriter)
        self.pushButtonAddLagrangianMesh.clicked.connect(self.slotAddLagrangianMesh)
        self.pushButtonDeleteLagrangianMesh.clicked.connect(self.slotDeleteLagrangianMesh)
        self.pushButtonAddAssociatedLagrangianWriter.clicked.connect(self.slotAddAssociatedLagrangianWriter)
        self.pushButtonDeleteAssociatedLagrangianWriter.clicked.connect(self.slotDeleteAssociatedLagrangianWriter)
        self.toolButtonAdd.clicked.connect(self.slotAddMonitoringPoint)
        self.toolButtonDelete.clicked.connect(self.slotDeleteMonitoringPoints)
        self.toolButtonDuplicate.clicked.connect(self.slotDuplicateMonitoringPoints)
        self.toolButtonImportCSV.clicked.connect(self.slotImportMonitoringPoints)

        # lineEdit
        self.lineEditNTLIST.textChanged[str].connect(self.slotListingFrequency)
        self.lineEditNTLAL.textChanged[str].connect(self.slotNTLAL)
        self.lineEditFrequency.textChanged[str].connect(self.slotWriterFrequency)
        self.lineEditFrequencyTime.textChanged[str].connect(self.slotWriterFrequencyTime)
        self.lineEditTimePlot.textChanged[str].connect(self.slotMonitoringPointFrequency)
        self.lineEditFRTimePlot.textChanged[str].connect(self.slotMonitoringPointFrequencyTime)

        # comboBox
        self.comboBoxOutput.activated[str].connect(self.slotOutputListing)
        self.comboBoxNTLAL.activated[str].connect(self.slotChoiceNTLAL)

        self.comboBoxTimeDependency.activated[str].connect(self.slotWriterTimeDependency)
        self.comboBoxFrequency.activated[str].connect(self.slotWriterFrequencyChoice)
        self.comboBoxFormatE.activated[str].connect(self.slotWriterOptions)
        self.comboBoxFormatH.activated[str].connect(self.slotWriterOptions)
        self.comboBoxPolygon.activated[str].connect(self.slotWriterOptions)
        self.comboBoxPolyhedra.activated[str].connect(self.slotWriterOptions)
        self.comboBoxTimePlot.activated[str].connect(self.slotMonitoringPoint)
        self.comboBoxProbeFmt.activated[str].connect(self.slotOutputProbeFmt)


        # Validators

        validatorNTLIST = IntValidator(self.lineEditNTLIST, min=1)
        validatorNTLAL = IntValidator(self.lineEditNTLAL)
        validatorFrequency = IntValidator(self.lineEditFrequency, min=1)
        validatorNTHIST = IntValidator(self.lineEditTimePlot, min=1)
        validatorFrequencyTime = DoubleValidator(self.lineEditFrequencyTime)
        validatorFRHIST = DoubleValidator(self.lineEditFRTimePlot)
        validatorRadius = DoubleValidator(self.lineEditProbesRadius, min=0.)
        validatorRadius.setExclusiveMin(True)

        self.lineEditNTLIST.setValidator(validatorNTLIST)
        self.lineEditNTLAL.setValidator(validatorNTLAL)
        self.lineEditFrequency.setValidator(validatorFrequency)
        self.lineEditTimePlot.setValidator(validatorNTHIST)
        self.lineEditFrequencyTime.setValidator(validatorFrequencyTime)
        self.lineEditFRTimePlot.setValidator(validatorFRHIST)
        self.lineEditProbesRadius.setValidator(validatorRadius)

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
        self.lineEditNTLIST.setText(str(ntlist))
        self.slotOutputListing(t)

        if self.case['run_type'] == 'standard':
            if self.lag_mdl.getLagrangianModel() != 'off':
                self.groupBoxListingParticles.show()
                period = self.mdl.getListingFrequencyLagrangian()
                if period == -1:
                    m = "None"
                elif period == 1:
                    m = "At each step"
                else:
                    m = "Frequency_l"
                self.modelNTLAL.setItem(str_model = m)
                t = self.modelNTLAL.dicoM2V[m]
                self.lineEditNTLAL.setText(str(period))
                self.slotChoiceNTLAL(t)

                # lagrangian mesh
                self.tabWidget.setTabEnabled(3, True)
            else:
                self.groupBoxListingParticles.hide()
                # lagrangian mesh
                self.tabWidget.setTabEnabled(3, False)
        else:
            self.groupBoxListingParticles.hide()
            # lagrangian mesh
            self.tabWidget.setTabEnabled(3, False)

        # Initialisation of the monitoring points files

        m = self.mdl.getMonitoringPointType()
        if m == 'Frequency_h_x' :
            frhist = self.mdl.getMonitoringPointFrequencyTime()
            self.lineEditFRTimePlot.setText(str(frhist))
        else :
            nthist = self.mdl.getMonitoringPointFrequency()
            self.lineEditTimePlot.setText(str(nthist))
        self.modelTimePlot.setItem(str_model=m)
        t = self.modelTimePlot.dicoM2V[m]
        self.slotMonitoringPoint(t)

        # Monitoring points initialisation

        if self.case['salome'] and not self.case['probes']:
            try:
                from SalomeActors import ProbeActors
                self.case['probes'] = ProbeActors()
                self.case['probes'].setTableView(self.tableViewPoints)
            except Exception:
                import traceback
                import sys
                exc_info = sys.exc_info()
                bt = traceback.format_exception(*exc_info)
                for l in bt:
                    print(l.rstrip())
                pass

        self.groupBoxProbesDisplay.setChecked(False)
        self.groupBoxProbesDisplay.setEnabled(False)

        if self.case['salome'] and self.case['probes']:
            self.case['probes'].removeActors()
            self.case['probes'].setTableView(self.tableViewPoints)
            self.groupBoxProbesDisplay.setChecked(self.case['probes'].getVisibility())
            self.groupBoxProbesDisplay.setEnabled(True)
            self.lineEditProbesRadius.setText(str(self.case['probes'].getRadius()))

            self.tableViewPoints.pressed[QModelIndex].connect(self.slotSelectedActors)
            self.tableViewPoints.entered[QModelIndex].connect(self.slotSelectedActors)
            self.groupBoxProbesDisplay.clicked[bool].connect(self.slotProbesDisplay)
            self.lineEditProbesRadius.textChanged.connect(self.slotProbesRadius)

        self.toolButtonDuplicate.setEnabled(False)
        if self.mdl.getNumberOfMonitoringPoints() > 0:
            self.toolButtonDuplicate.setEnabled(True)

        for n in range(self.mdl.getNumberOfMonitoringPoints()):
            name = str(n+1)
            X, Y, Z = self.mdl.getMonitoringPointCoordinates(name)
            self.__insertMonitoringPoint(name, X, Y, Z)
            if self.case['salome']:
                self.__salomeHandlerAddMonitoringPoint(name, X, Y, Z)

        # Writer initialisation

        self.groupBoxFrequency.hide()
        self.groupBoxTimeDependency.hide()
        self.groupBoxOptions.hide()

        # Mesh initialisation

        self.groupBoxVariable.hide()
        self.groupBoxAssociatedWriter.hide()

        # Mesh initialisation
        self.groupBoxLagrangianVariable.hide()
        self.groupBoxAssociatedLagrangianWriter.hide()

        # values of probes format
        fmt = self.mdl.getMonitoringPointFormat()
        self.modelProbeFmt.setItem(str_model=fmt)

        # tabWidget active
        self.tabWidget.setCurrentIndex(self.case['current_tab'])

        self.case.undoStartGlobal()


    @pyqtSlot(str)
    def slotOutputListing(self, text):
        """
        INPUT choice of the output listing
        """
        listing = self.modelOutput.dicoV2M[str(text)]
        log.debug("slotOutputListing-> listing = %s" % listing)

        if listing == "None":
            ntlist = -1
            self.mdl.setListingFrequency(ntlist)
            self.lineEditNTLIST.setText(str(ntlist))
            self.lineEditNTLIST.setDisabled(True)

        elif listing == "At each step":
            ntlist = 1
            self.lineEditNTLIST.setText(str(ntlist))
            self.lineEditNTLIST.setDisabled(True)

        elif listing == "Frequency_l":
            self.lineEditNTLIST.setEnabled(True)
            ntlist = from_qvariant(self.lineEditNTLIST.text(), int)
            if ntlist < 1:
                ntlist = 1
                self.lineEditNTLIST.setText(str(ntlist))


    @pyqtSlot(str)
    def slotChoiceNTLAL(self, text):
        """
        INPUT choice of the output listing for lagrangian variables
        """
        listing = self.modelNTLAL.dicoV2M[str(text)]
        log.debug("slotChoiceNTLAL-> listing = %s" % listing)

        if listing == "None":
            ntlist = -1
            self.mdl.setListingFrequencyLagrangian(ntlist)
            self.lineEditNTLAL.setText(str(ntlist))
            self.lineEditNTLAL.setDisabled(True)

        elif listing == "At each step":
            ntlist = 1
            self.mdl.setListingFrequencyLagrangian(ntlist)
            self.lineEditNTLAL.setText(str(ntlist))
            self.lineEditNTLAL.setDisabled(True)

        elif listing == "Frequency_l":
            self.lineEditNTLAL.setEnabled(True)
            ntlist = from_qvariant(self.lineEditNTLAL.text(), int)
            if ntlist < 1:
                ntlist = 1
                self.mdl.setListingFrequencyLagrangian(ntlist)
                self.lineEditNTLAL.setText(str(ntlist))


    @pyqtSlot(str)
    def slotListingFrequency(self, text):
        """
        Input the frequency of the listing output
        """
        if self.lineEditNTLIST.validator().state == QValidator.Acceptable:
            n = from_qvariant(text, int)
            log.debug("slotListingFrequency-> NTLIST = %s" % n)
            self.mdl.setListingFrequency(n)


    @pyqtSlot(str)
    def slotNTLAL(self, text):
        """
        Input the frequency of the listing output for lagrangian variables
        """
        if self.lineEditNTLAL.validator().state == QValidator.Acceptable:
            n = from_qvariant(text, int)
            log.debug("slotNTLAL-> NTLIST = %s" % n)
            self.mdl.setListingFrequencyLagrangian(n)


    def __insertWriter(self, name, writer_id, format, directory):
        """
        Add a new 'item' into the Hlist.
        """
        self.modelWriter.newData(name, writer_id, format, directory)


    @pyqtSlot()
    def slotAddWriter(self):
        """
        Add one monitoring point with these coordinates in the list in the Hlist
        The number of the monitoring points is added to the preceding one
        """
        writer_id = self.mdl.addWriter()
        self.__insertWriter(self.mdl.getWriterLabel(writer_id),
                            writer_id,
                            self.mdl.getWriterFormat(writer_id),
                            self.mdl.getWriterDirectory(writer_id))

    @pyqtSlot()
    def slotDeleteWriter(self):
        """
        Just delete the current selected entries from the Hlist and
        of course from the XML file.
        """
        lst = []
        selectionModel = self.tableViewWriter.selectionModel()
        for index in selectionModel.selectedRows():
            w = self.modelWriter.getItem(index.row())['id']
            if int(w) < 0:
                title = self.tr("Warning")
                msg   = self.tr("You can't delete a default writer.")
                QMessageBox.information(self, title, msg)
                return
            lst.append(str(w))

        self.mdl.deleteWriter(lst)

        self.modelWriter.deleteAllData()
        list_writer = []
        for writer in self.mdl.getWriterIdList():
            if int(writer) > 0:
                list_writer.append(writer)
        for writer in self.mdl.getWriterIdList():
            if int(writer) < 0:
                label = self.mdl.getWriterLabel(writer)
                format = self.mdl.getWriterFormat(writer)
                directory = self.mdl.getWriterDirectory(writer)
                self.__insertWriter(label, writer, format, directory)
        new_id = 0
        for writer in list_writer:
            new_id = new_id + 1
            label = self.mdl.getWriterLabel(writer)
            format = self.mdl.getWriterFormat(writer)
            directory = self.mdl.getWriterDirectory(writer)
            self.__insertWriter(label, str(new_id), format, directory)
        self.tableViewMesh.clearSelection()
        self.groupBoxVariable.hide()
        self.groupBoxAssociatedWriter.hide()
        self.groupBoxLagrangianVariable.hide()
        self.groupBoxAssociatedLagrangianWriter.hide()

    @pyqtSlot(QModelIndex,QModelIndex)
    def slotDataChanged(self, topLeft, bottomRight):
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


    def showAssociatedWriterTable(self):
        cindex = self.tableViewMesh.currentIndex()
        if cindex != (-1,-1):
            row = cindex.row()
            mesh_id = self.modelMesh.getItem(row)['id']

            self.modelAssociatedWriter = StandardItemModelAssociatedWriter(self, self.mdl, mesh_id)
            if QT_API == "PYQT4":
                self.tableViewAssociatedWriter.horizontalHeader().setResizeMode(QHeaderView.Stretch)
            elif QT_API == "PYQT5":
                self.tableViewAssociatedWriter.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)

            delegate_associated_writer = AssociatedWriterDelegate(self.tableViewAssociatedWriter, 0)
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


    def showAssociatedLagrangianWriterTable(self):
        cindex = self.tableViewLagrangianMesh.currentIndex()
        if cindex != (-1,-1):
            row = cindex.row()
            mesh_id = self.modelLagrangianMesh.getItem(row)['id']

            self.modelAssociatedLagrangianWriter = StandardItemModelAssociatedWriter(self, self.mdl, mesh_id)
            if QT_API == "PYQT4":
                self.tableViewAssociatedLagrangianWriter.horizontalHeader().setResizeMode(QHeaderView.Stretch)
            elif QT_API == "PYQT5":
                self.tableViewAssociatedLagrangianWriter.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)

            delegate_associated_writer = AssociatedWriterDelegate(self.tableViewAssociatedLagrangianWriter, 1)
            self.tableViewAssociatedLagrangianWriter.setItemDelegateForColumn(0, delegate_associated_writer)

            self.tableViewAssociatedLagrangianWriter.reset()
            self.modelAssociatedLagrangianWriter = StandardItemModelAssociatedWriter(self, self.mdl, mesh_id)
            self.tableViewAssociatedLagrangianWriter.setModel(self.modelAssociatedLagrangianWriter)
            self.modelAssociatedLagrangianWriter.deleteAllData()
            writer_row = 0
            for n in self.mdl.getAssociatedWriterIdList(mesh_id):
                label = self.mdl.getWriterLabel(n)
                self.__insertAssociatedLagrangianWriter(label)
                writer_row = writer_row +1


    @pyqtSlot("QModelIndex")
    def slotSelectWriter(self, index):
        cindex = self.tableViewWriter.currentIndex()
        if cindex != (-1,-1):
            row = cindex.row()
            writer_id = self.modelWriter.getItem(row)['id']
            self.groupBoxFrequency.show()
            self.groupBoxTimeDependency.show()
            self.groupBoxOptions.show()

            if writer_id == "-3" or writer_id == "-4":
                self.comboBoxTimeDependency.setEnabled(False)
            else:
                self.comboBoxTimeDependency.setEnabled(True)

            frequency_choice = self.mdl.getWriterFrequencyChoice(writer_id)
            self.modelFrequency.setItem(str_model=frequency_choice)

            if frequency_choice == "none":
                self.lineEditFrequency.hide()
                self.lineEditFrequencyTime.hide()
                self.pushButtonFrequency.hide()

            if frequency_choice == "time_step":
                self.lineEditFrequency.show()
                self.lineEditFrequency.setEnabled(True)
                ntchr = int(self.mdl.getWriterFrequency(writer_id))
                if ntchr < 1:
                    ntchr = 1
                    self.mdl.setWriterFrequency(writer_id, ntchr)
                self.lineEditFrequency.setText(str(ntchr))
                self.lineEditFrequencyTime.hide()
                self.pushButtonFrequency.hide()

            if frequency_choice == "time_value":
                self.lineEditFrequency.hide()
                self.lineEditFrequencyTime.show()
                frchr = float(self.mdl.getWriterFrequency(writer_id))
                self.lineEditFrequencyTime.setText(str(frchr))
                self.pushButtonFrequency.hide()

            if frequency_choice == "formula":
                self.lineEditFrequency.hide()
                self.lineEditFrequencyTime.hide()
                self.pushButtonFrequency.show()
                self.pushButtonFrequency.setEnabled(True)
                exp = self.mdl.getWriterFrequency(writer_id)
                if exp:
                    self.pushButtonFrequency.setStyleSheet("background-color: green")
                    self.pushButtonFrequency.setToolTip(exp)
                else:
                    self.pushButtonFrequency.setStyleSheet("background-color: red")

            if self.mdl.getWriterOutputStartStatus(writer_id) == 'on':
                self.checkBoxOutputStart.setChecked(True)
            else:
                self.checkBoxOutputStart.setChecked(False)

            if self.mdl.getWriterOutputEndStatus(writer_id) == 'on':
                self.checkBoxOutputEnd.setChecked(True)
            else:
                self.checkBoxOutputEnd.setChecked(False)

            time_dependency = self.mdl.getWriterTimeDependency(writer_id)
            self.modelTimeDependency.setItem(str_model=time_dependency)
            options = self.mdl.getWriterOptions(writer_id)
            self.__updateOptionsFormat(options, row)


    @pyqtSlot(str)
    def slotWriterFrequencyChoice(self, text):
        """
        INPUT choice of the output frequency for a writer
        """
        cindex = self.tableViewWriter.currentIndex()
        if cindex != (-1, -1):
            row = cindex.row()
            writer_id = self.modelWriter.getItem(row)['id']
            chrono = self.modelFrequency.dicoV2M[str(text)]
            log.debug("slotOutputPostpro-> chrono = %s" % chrono)
            self.mdl.setWriterFrequencyChoice(writer_id, chrono)

            if chrono == "none":
                self.lineEditFrequency.hide()
                self.lineEditFrequencyTime.hide()
                self.pushButtonFrequency.hide()

            elif chrono == "time_step":
                self.lineEditFrequency.show()
                self.lineEditFrequency.setEnabled(True)
                self.pushButtonFrequency.setEnabled(False)
                ntchr = int(self.mdl.getWriterFrequency(writer_id))
                if ntchr < 1:
                    ntchr = 1
                    self.mdl.setWriterFrequency(writer_id, ntchr)
                self.lineEditFrequency.setText(str(ntchr))
                self.lineEditFrequencyTime.hide()
                self.pushButtonFrequency.hide()

            elif chrono == "time_value":
                self.lineEditFrequency.hide()
                self.lineEditFrequencyTime.show()
                self.pushButtonFrequency.setEnabled(False)
                frchr = self.mdl.getWriterFrequency(writer_id)
                self.lineEditFrequencyTime.setText(str(frchr))
                self.pushButtonFrequency.hide()

            elif chrono == "formula":
                self.lineEditFrequency.hide()
                self.lineEditFrequencyTime.hide()
                self.pushButtonFrequency.setEnabled(True)
                exp = self.mdl.getWriterFrequency(writer_id)
                if exp:
                    self.pushButtonFrequency.setStyleSheet("background-color: green")
                    self.pushButtonFrequency.setToolTip(exp)
                else:
                    self.pushButtonFrequency.setStyleSheet("background-color: red")
                self.pushButtonFrequency.show()



    @pyqtSlot(str)
    def slotWriterFrequency(self, text):
        """
        Input the frequency of the post-processing output
        """
        cindex = self.tableViewWriter.currentIndex()
        if cindex != (-1,-1):
            row = cindex.row()
            writer_id = self.modelWriter.getItem(row)['id']
            self.lineEditFrequency.setEnabled(True)
            # Robustness to avoid crash when text is empty:
            if self.lineEditFrequency.text() != '':
                n = from_qvariant(self.lineEditFrequency.text(), int)
            else:
                n = 1
            if self.lineEditFrequency.validator().state == QValidator.Acceptable:
                log.debug("slotPostproFrequency = %s" % n)
                self.mdl.setWriterFrequency(writer_id, str(n))


    @pyqtSlot(str)
    def slotWriterFrequencyTime(self, text):
        """
        Input the frequency of the post-processing output
        """
        cindex = self.tableViewWriter.currentIndex()
        if cindex != (-1,-1):
            row = cindex.row()
            writer_id = self.modelWriter.getItem(row)['id']
            if self.lineEditFrequencyTime.validator().state == QValidator.Acceptable:
                n = from_qvariant(text, float)
                log.debug("slotPostproFrequencyTime-> FRCHR = %s" % n)
                self.mdl.setWriterFrequency(writer_id, str(n))


    @pyqtSlot()
    def slotWriterFrequencyFormula(self):
        """
        """
        cindex = self.tableViewWriter.currentIndex()
        if cindex != (-1,-1):
            row = cindex.row()
            writer_id = self.modelWriter.getItem(row)['id']
            exp = self.mdl.getWriterFrequency(writer_id)
            if not exp:
                exp = """iactive = 1;\n"""
            exa = """#example:"""
            req = [('iactive', 'at a time step the writer is active or not')]
            sym = [('t', 'current time'),
                   ('niter', 'current time step')]

            for (nme, val) in self.notebook.getNotebookList():
                sym.append((nme, 'value (notebook) = ' + str(val)))

            dialog = QMeiEditorView(self,
                                    check_syntax = self.case['package'].get_check_syntax(),
                                    expression = exp,
                                    required   = req,
                                    symbols    = sym,
                                    examples   = exa)
            if dialog.exec_():
                result = str(dialog.get_result())
                log.debug("slotWriterFrequencyFormula -> %s" % result)
                self.mdl.setWriterFrequency(writer_id, result)
                self.pushButtonFrequency.setStyleSheet("background-color: green")
                self.pushButtonFrequency.setToolTip(result)


    @pyqtSlot(str)
    def slotWriterTimeDependency(self, text):
        """
        Input type of post-processing for mesh
        """
        cindex = self.tableViewWriter.currentIndex()
        if cindex != (-1,-1):
            row = cindex.row()
            writer_id = self.modelWriter.getItem(row)['id']
            self.mdl.setWriterTimeDependency(writer_id,
                                             self.modelTimeDependency.dicoV2M[str(text)])


    @pyqtSlot()
    def slotWriterOutputStart(self):
        """
        Input output start flag
        """
        cindex = self.tableViewWriter.currentIndex()
        if cindex != (-1,-1):
            row = cindex.row()
            writer_id = self.modelWriter.getItem(row)['id']
            st = 'on'
            if not self.checkBoxOutputStart.isChecked():
              st = 'off'
            self.mdl.setWriterOutputStartStatus(writer_id, st)


    @pyqtSlot()
    def slotWriterOutputEnd(self):
        """
        Input output end flag
        """
        cindex = self.tableViewWriter.currentIndex()
        if cindex != (-1,-1):
            row = cindex.row()
            writer_id = self.modelWriter.getItem(row)['id']
            st = 'on'
            if not self.checkBoxOutputEnd.isChecked():
              st = 'off'
            self.mdl.setWriterOutputEndStatus(writer_id, st)


    def __WriterOptionsPrelude(self):
        """
        Common handling for slots handling writer options
        """
        row = -1
        writer_id = None
        options = []
        cindex = self.tableViewWriter.currentIndex()
        if cindex != (-1,-1):
            row = cindex.row()
            writer_id = self.modelWriter.getItem(row)['id']
            l = re.split(r'[;,\s]\s*',
                         self.mdl.getWriterOptions(writer_id))
            for s in l:
                if s:
                    options.append(s)

        return row, writer_id, options


    @pyqtSlot()
    def slotWriterSeparateMeshes(self):
        """
        Writer separate meshes option
        """
        row, writer_id, options = self.__WriterOptionsPrelude()

        if not writer_id:  # should not occur
            return

        if self.checkBoxSeparateMeshes.isChecked():
            if 'separate_meshes' not in options:
                options.append('separate_meshes')
        else:
            if 'separate_meshes' in options:
                options.remove('separate_meshes')

        l = ', '.join(options)
        log.debug("slotWriterSeparateMeshes")
        self.mdl.setWriterOptions(writer_id, l)


    @pyqtSlot()
    def slotWriterParallelIO(self):
        """
        Writer separate meshes option
        """
        row, writer_id, options = self.__WriterOptionsPrelude()

        if not writer_id:  # should not occur
            return

        if self.checkBoxAllowParallelIO.isChecked():
            if 'serial_io' in options:
                options.remove('serial_io')
        else:
            if 'serial_io' not in options:
                options.append('serial_io')

        l = ', '.join(options)
        log.debug("slotWriterParallelIO")
        self.mdl.setWriterOptions(writer_id, l)


    @pyqtSlot(str)
    def slotWriterOptions(self):
        """
        Create line for command of format's options
        """
        row, writer_id, options = self.__WriterOptionsPrelude()

        if not writer_id:  # should not occur
            return

        writer_id = self.modelWriter.getItem(row)['id']
        line = []

        # remove previous settings matching this slot
        for m in (self.modelFormatE, self.modelFormatH,
                  self.modelPolygon, self.modelPolyhedra):
            for s in m.items:
                if s in options:
                    options.remove(s)

        # now append options
        format = self.modelWriter.getItem(row)['format']
        if format == 'ensight':
            opt_format = self.modelFormatE.dicoV2M[str(self.comboBoxFormatE.currentText())]
            if opt_format != 'binary':
                options.append(opt_format)
        elif format == 'histogram':
            opt_format = self.modelFormatH.dicoV2M[str(self.comboBoxFormatH.currentText())]
            if opt_format != 'text':
                options.append(opt_format)

        opt_polygon = self.modelPolygon.dicoV2M[str(self.comboBoxPolygon.currentText())]
        opt_polyhed = self.modelPolyhedra.dicoV2M[str(self.comboBoxPolyhedra.currentText())]
        if opt_polygon != 'display':
            options.append(opt_polygon)
        if opt_polyhed != 'display':
            options.append(opt_polyhed)

        l = ', '.join(options)
        log.debug("slotOutputOptions-> %s" % l)
        self.mdl.setWriterOptions(writer_id, l)


    def __updateOptionsFormat(self, options, row):
        """
        Update line for command of format's options at each modification of
        post processing format
        """
        opts = re.split(r'[;,\s]\s*', options)

        format = self.modelWriter.getItem(row)['format']
        log.debug("__updateOptionsFormat-> format = %s" % format)
        log.debug("__updateOptionsFormat-> options = %s" % options)

        # update widgets from the options list

        parallel_io = True
        for opt in opts:
            if format == 'ensight' and opt in ['binary', 'big_endian', 'text']:
                self.modelFormatE.setItem(str_model=opt)
            elif format == 'histogram' and opt in ['text', 'tex', 'png']:
                self.modelFormatH.setItem(str_model=opt)
            elif opt == 'discard_polygons' or opt == 'divide_polygons':
                self.modelPolygon.setItem(str_model=opt)
            elif opt == 'discard_polyhedra' or opt == 'divide_polyhedra':
                self.modelPolyhedra.setItem(str_model=opt)
            elif opt == 'separate_meshes':
                self.checkBoxSeparateMeshes.setChecked(True)
            elif opt == 'serial_io':
                parallel_io = False

        # default

        if 'discard_polygons' not in opts and 'divide_polygons' not in opts:
            self.modelPolygon.setItem(str_model="display")
        if 'discard_polyhedra' not in opts and 'divide_polyhedra' not in opts:
            self.modelPolyhedra.setItem(str_model="display")
        if 'separate_meshes' not in opts:
            self.checkBoxSeparateMeshes.setChecked(False)

        # enable and disable options related to the format

        if format in ["ensight", "med", "cgns"]:
            self.modelPolygon.enableItem(str_model='discard_polygons')
            self.modelPolygon.enableItem(str_model='divide_polygons')
            self.comboBoxPolygon.show()
            self.labelPolygon.show()
            if format == "cgns":
                self.modelPolyhedra.setItem(str_model='divide_polyhedra')
                self.modelPolyhedra.disableItem(str_model='display')
            else:
                self.modelPolyhedra.enableItem(str_model='divide_polyhedra')
            self.modelPolyhedra.enableItem(str_model='discard_polyhedra')
            if format == "med":
                self.checkBoxAllowParallelIO.show()
                self.checkBoxAllowParallelIO.setChecked(parallel_io)
            else:
                self.checkBoxAllowParallelIO.hide()
            self.comboBoxPolyhedra.show()
            self.labelPolyhedra.show()
        else:
            self.checkBoxAllowParallelIO.hide()
            self.comboBoxPolygon.hide()
            self.labelPolygon.hide()
            self.comboBoxPolyhedra.hide()
            self.labelPolyhedra.hide()

        if format == "ensight":
            if 'binary' not in opts and 'big_endian' not in opts and 'text' not in opts:
                self.modelFormatE.setItem(str_model='binary')
            self.modelFormatE.show()
            self.labelFormatE.show()
        else:
            self.modelFormatE.hide()
            self.labelFormatE.hide()

        if format == "histogram":
            if 'text' not in opts and 'tex' not in opts and 'png' not in opts:
                self.modelFormatH.setItem(str_model='text')
            self.modelFormatH.show()
            self.labelFormatH.show()
        else:
            self.modelFormatH.hide()
            self.labelFormatH.hide()


    def __insertMesh(self, name, mesh_id, mesh_type, selection):
        """
        Add a new 'item' into the Hlist.
        """
        self.modelMesh.newData(name, mesh_id, mesh_type, selection)


    def __insertLagrangianMesh(self, name, mesh_id, mesh_type, density, selection):
        """
        Add a new 'item' into the Hlist.
        """
        self.modelLagrangianMesh.newData(name, mesh_id, mesh_type, density, selection)


    @pyqtSlot()
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


    @pyqtSlot()
    def slotDeleteMesh(self):
        """
        Just delete the current selected entries from the Hlist and
        of course from the XML file.
        """
        lst = []
        selectionModel = self.tableViewMesh.selectionModel()
        for index in selectionModel.selectedRows():
            mesh_id = self.modelMesh.getItem(index.row())['id']
            if int(mesh_id) < 0:
                title = self.tr("Warning")
                msg   = self.tr("You can't delete a default mesh\n"
                                "(but you may disassociate it from all writers).")
                QMessageBox.information(self, title, msg)
                return
            lst.append(str(mesh_id))

        self.mdl.deleteMesh(lst)

        self.modelMesh.deleteAllData()
        self.modelLagrangianMesh.deleteAllData()
        list_mesh = []
        for mesh in self.mdl.getMeshIdList():
            if int(mesh) > 0:
                list_mesh.append(mesh)
        new_id = 0
        for mesh in self.mdl.getMeshIdList():
            if int(mesh) < 0:
                label = self.mdl.getMeshLabel(mesh)
                mesh_type = self.mdl.getMeshType(mesh)
                location = self.mdl.getMeshLocation(mesh)
                if mesh_type != "particles":
                    self.__insertMesh(label, mesh, mesh_type, location)
                else:
                    density = self.mdl.getMeshDensity(mesh)
                    self.__insertLagrangianMesh(label, mesh, mesh_type, density, location)
        for mesh in list_mesh:
            new_id = new_id + 1
            label = self.mdl.getMeshLabel(mesh)
            mesh_type = self.mdl.getMeshType(mesh)
            location = self.mdl.getMeshLocation(mesh)
            if mesh_type != "particles":
                self.__insertMesh(label, str(new_id), mesh_type, location)
            else:
                density = self.mdl.getMeshDensity(mesh)
                self.__insertLagrangianMesh(label, str(new_id), mesh_type, density, location)


    @pyqtSlot()
    def slotAddLagrangianMesh(self):
        """
        Add one monitoring point with these coordinates in the list in the Hlist
        The number of the monitoring point is added at the precedent one
        """
        mesh_id = self.mdl.addLagrangianMesh()
        self.__insertLagrangianMesh(self.mdl.getMeshLabel(mesh_id),
                                    mesh_id,
                                    self.mdl.getLagrangianMeshType(mesh_id),
                                    self.mdl.getMeshDensity(mesh_id),
                                    self.mdl.getMeshLocation(mesh_id))


    @pyqtSlot()
    def slotDeleteLagrangianMesh(self):
        """
        Just delete the current selected entries from the Hlist and
        of course from the XML file.
        """
        lst = []
        selectionModel = self.tableViewLagrangianMesh.selectionModel()
        for index in selectionModel.selectedRows():
            mesh_id = self.modelLagrangianMesh.getItem(index.row())['id']
            if int(mesh_id) < 0:
                title = self.tr("Warning")
                msg   = self.tr("You can't delete a default mesh\n"
                                "(but you may disassociate it from all writers).")
                QMessageBox.information(self, title, msg)
                return
            lst.append(str(mesh_id))

        self.mdl.deleteMesh(lst)

        self.modelMesh.deleteAllData()
        self.modelLagrangianMesh.deleteAllData()
        list_mesh = []
        for mesh in self.mdl.getMeshIdList():
            if int(mesh) > 0:
                list_mesh.append(mesh)
        new_id = 0
        for mesh in self.mdl.getMeshIdList():
            if int(mesh) < 0:
                label = self.mdl.getMeshLabel(mesh)
                mesh_type = self.mdl.getMeshType(mesh)
                location = self.mdl.getMeshLocation(mesh)
                if mesh_type != "particles":
                    self.__insertMesh(label, mesh, mesh_type, location)
                else:
                    density = self.mdl.getMeshDensity(mesh)
                    self.__insertLagrangianMesh(label, mesh, mesh_type, density, location)
        for mesh in list_mesh:
            new_id = new_id + 1
            label = self.mdl.getMeshLabel(mesh)
            mesh_type = self.mdl.getMeshType(mesh)
            location = self.mdl.getMeshLocation(mesh)
            if mesh_type != "particles":
                self.__insertMesh(label, str(new_id), mesh_type, location)
            else:
                density = self.mdl.getMeshDensity(mesh)
                self.__insertLagrangianMesh(label, str(new_id), mesh_type, density, location)


    @pyqtSlot("QModelIndex")
    def slotSelectMesh(self, index):
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


    @pyqtSlot("QModelIndex")
    def slotSelectLagrangianMesh(self, index):
        cindex = self.tableViewLagrangianMesh.currentIndex()
        if cindex != (-1,-1):
            row = cindex.row()
            mesh_id = self.modelLagrangianMesh.getItem(row)['id']
            self.groupBoxLagrangianVariable.show()
            self.groupBoxAssociatedLagrangianWriter.show()
            if int(mesh_id) < 0:
                self.checkBoxAllLagrangianVariables.setEnabled(False)
                self.checkBoxAllLagrangianVariables.setChecked(True)
                self.mdl.setMeshAllVariablesStatus(mesh_id,"on")
            else:
                self.checkBoxAllLagrangianVariables.setEnabled(True)
                all_variables = self.mdl.getMeshAllVariablesStatus(mesh_id)
                if all_variables == 'on':
                    self.checkBoxAllLagrangianVariables.setChecked(True)
                else :
                    self.checkBoxAllLagrangianVariables.setChecked(False)
            self.showAssociatedLagrangianWriterTable()


    @pyqtSlot()
    def slotAllVariables(self):
        """
        Toggle all variables output
        """
        cindex = self.tableViewMesh.currentIndex()
        if cindex != (-1,-1):
            row = cindex.row()
            mesh_id = self.modelMesh.getItem(row)['id']
            if self.checkBoxAllVariables.isChecked():
                self.mdl.setMeshAllVariablesStatus(mesh_id, "on")
            else:
                self.mdl.setMeshAllVariablesStatus(mesh_id, "off")


    @pyqtSlot()
    def slotAllLagrangianVariables(self):
        """
        Toggle all lagrangian variables output
        """
        cindex = self.tableViewLagrangianMesh.currentIndex()
        if cindex != (-1,-1):
            row = cindex.row()
            mesh_id = self.modelLagrangianMesh.getItem(row)['id']
            if self.checkBoxAllLagrangianVariables.isChecked():
                self.mdl.setMeshAllVariablesStatus(mesh_id, "on")
            else:
                self.mdl.setMeshAllVariablesStatus(mesh_id, "off")


    def __insertAssociatedWriter(self, name):
        """
        Add a new 'item' into the Hlist.
        """
        self.modelAssociatedWriter.newItem(name)


    def __insertAssociatedLagrangianWriter(self, name):
        """
        Add a new 'item' into the Hlist.
        """
        self.modelAssociatedLagrangianWriter.newItem(name)


    @pyqtSlot()
    def slotAddAssociatedWriter(self):
        """
        Add one monitoring point with these coordinates in the list in the Hlist
        The number of the monitoring point is added at the precedent one
        """
        cindex = self.tableViewMesh.currentIndex()
        if cindex != (-1,-1):
            row = cindex.row()
            mesh_id = self.modelMesh.getItem(row)['id']
            lagrangian = 0
            associated_writer_id = self.mdl.addAssociatedWriter(mesh_id, lagrangian)
            if associated_writer_id == None:
                title = self.tr("Warning")
                msg   = self.tr("Please create another writer\n"\
                                "before adding a new associated writer.")
                QMessageBox.information(self, title, msg)
                return
            self.__insertAssociatedWriter(self.mdl.getWriterLabel(associated_writer_id))


    @pyqtSlot()
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


    @pyqtSlot()
    def slotAddAssociatedLagrangianWriter(self):
        """
        Add one monitoring point with these coordinates in the list in the Hlist
        The number of the monitoring point is added at the precedent one
        """
        cindex = self.tableViewLagrangianMesh.currentIndex()
        if cindex != (-1,-1):
            row = cindex.row()
            mesh_id = self.modelLagrangianMesh.getItem(row)['id']
            lagrangian = 1
            associated_writer_id = self.mdl.addAssociatedWriter(mesh_id, lagrangian)
            if associated_writer_id == None:
                title = self.tr("Warning")
                msg   = self.tr("Please create another writer\n"\
                                "before adding a new associated writer.")
                QMessageBox.information(self, title, msg)
                return
            self.__insertAssociatedLagrangianWriter(self.mdl.getWriterLabel(associated_writer_id))


    @pyqtSlot()
    def slotDeleteAssociatedLagrangianWriter(self):
        """
        Just delete the current selected entries from the Hlist and
        of course from the XML file.
        """
        cindex = self.tableViewLagrangianMesh.currentIndex()
        if cindex != (-1,-1):
            row = cindex.row()
            mesh_id = self.modelLagrangianMesh.getItem(row)['id']
            selectionModel = self.tableViewAssociatedLagrangianWriter.selectionModel()
            for index in selectionModel.selectedRows():
                writer_label = self.modelAssociatedLagrangianWriter.getItem(index.row())
                writer_id = self.mdl.getWriterIdFromLabel(writer_label)
                self.mdl.deleteAssociatedWriter(mesh_id, writer_id)

                self.modelAssociatedLagrangianWriter.deleteAllData()
                list_associated_writer = []
                for associated_writer in self.mdl.getAssociatedWriterIdList(mesh_id):
                    list_associated_writer.append(associated_writer)
                for associated_writer in list_associated_writer:
                    label = self.mdl.getWriterLabel(associated_writer)
                    self.__insertAssociatedLagrangianWriter(label)


    @pyqtSlot(str)
    def slotOutputProbeFmt(self, text):
        """
        INPUT choice of the output for the probes (.dat, .csv)
        """
        fmt = self.modelProbeFmt.dicoV2M[str(text)]
        log.debug("slotOutputProbeFmt-> fmt = %s" % fmt)
        self.mdl.setMonitoringPointFormat(fmt)


    @pyqtSlot(str)
    def slotMonitoringPoint(self, text):
        """
        Input choice of the output of monitoring points files
        """
        histo = self.modelTimePlot.dicoV2M[str(text)]
        log.debug("slotMonitoringPoint-> histo = %s" % histo)
        self.mdl.setMonitoringPointType(histo)

        if histo == "None":
            nthist = -1
            self.mdl.setMonitoringPointFrequency(nthist)
            self.lineEditTimePlot.hide()
            self.lineEditFRTimePlot.hide()
            self.comboBoxProbeFmt.hide()
            self.label.hide()
        else:
            self.comboBoxProbeFmt.show()
            self.label.show()

        if histo == "At each step":
            nthist = 1
            self.mdl.setMonitoringPointFrequency(nthist)
            self.lineEditTimePlot.hide()
            self.lineEditFRTimePlot.hide()

        if histo == "Frequency_h":
            self.lineEditTimePlot.show()
            self.lineEditTimePlot.setEnabled(True)
            nthist = self.mdl.getMonitoringPointFrequency()
            if nthist < 1:
                nthist = 1
                self.mdl.setMonitoringPointFrequency(nthist)
            self.lineEditTimePlot.setText(str(nthist))
            self.lineEditFRTimePlot.hide()

        if histo == "Frequency_h_x":
            self.lineEditTimePlot.hide()
            self.lineEditFRTimePlot.show()
            frlist = self.mdl.getMonitoringPointFrequencyTime()
            self.lineEditFRTimePlot.setText(str(frlist))


    @pyqtSlot(str)
    def slotMonitoringPointFrequencyTime(self, text):
        """
        Input the frequency of the monitoring point output
        """
        if self.lineEditFRTimePlot.validator().state == QValidator.Acceptable:
            n = from_qvariant(text, float)
            log.debug("slotMonitoringPointFrequencyTime-> FRHIST = %s" % n)
            self.mdl.setMonitoringPointFrequencyTime(n)


    @pyqtSlot(str)
    def slotMonitoringPointFrequency(self, text):
        """
        Input the frequency of the monitoring point output
        """
        if self.lineEditTimePlot.validator().state == QValidator.Acceptable:
            n = from_qvariant(text, int)
            log.debug("slotMonitoringPointFrequency-> NTHIST = %s" % n)
            self.mdl.setMonitoringPointFrequency(n)


    def __insertMonitoringPoint(self, num, X, Y, Z):
        """
        Add a new 'item' into the Hlist.
        """
        self.modelMonitoring.insertData(num, X, Y, Z)


    @pyqtSlot()
    def slotAddMonitoringPoint(self):
        """
        Add one monitoring point with these coordinates in the list in the Hlist
        The number of the monitoring point is added at the precedent one
        """
        self.mdl.addMonitoringPoint(x=0.0, y=0.0, z=0.0)
        n = self.mdl.getNumberOfMonitoringPoints()
        self.__insertMonitoringPoint(n, str('0'), str('0'), str('0'))

        self.toolButtonDuplicate.setEnabled(True)

        if self.case['salome']:
            self.__salomeHandlerAddMonitoringPoint(n, 0., 0., 0.)


    @pyqtSlot()
    def slotDeleteMonitoringPoints(self):
        """
        Just delete the current selected entries from the Hlist and
        of course from the XML file.
        If salome, delete from the VTK view.
        """
        l1 = []
        l2 = []
        selectionModel = self.tableViewPoints.selectionModel()
        for index in selectionModel.selectedRows():
            name = index.row() + 1
            l1.append(name)
            l2.append(name)

        log.debug("slotDeleteMonitoringPoints -> %s" % (l1,))

        self.mdl.deleteMonitoringPoints(l1)

        self.modelMonitoring.deleteAllData()
        for n in range(self.mdl.getNumberOfMonitoringPoints()):
            name = str(n + 1)
            X, Y, Z = self.mdl.getMonitoringPointCoordinates(name)
            self.__insertMonitoringPoint(name, X, Y, Z)

        if self.case['salome']:
            self.__salomeHandlerDeleteMonitoringPoint(l2)

        if self.mdl.getNumberOfMonitoringPoints() == 0:
            self.toolButtonDuplicate.setEnabled(False)


    @pyqtSlot()
    def slotDuplicateMonitoringPoints(self):
        """
        Duplicate monitoring points selected with these coordinates in the list in the Hlist
        """
        log.debug("slotDuplicateMonitoringPoints")

        l1 = []
        l2 = []
        selectionModel = self.tableViewPoints.selectionModel()

        probe_number = self.mdl.getNumberOfMonitoringPoints()

        idx = 1

        for index in selectionModel.selectedRows():
            name = str(index.row() + 1)
            X, Y, Z = self.mdl.getMonitoringPointCoordinates(name)
            new_name = str(probe_number + idx)
            self.mdl.addMonitoringPoint(x=X, y=Y, z=Z)
            self.__insertMonitoringPoint(new_name, X, Y, Z)

            if self.case['salome']:
                self.__salomeHandlerAddMonitoringPoint(probe_number + idx, X, Y, Z)

            idx = idx + 1


    @pyqtSlot()
    def slotImportMonitoringPoints(self):
        """
        select a csv file to add probes
        """
        log.debug("slotImportMonitoringPoints")

        data = self.case['data_path']
        title = self.tr("Probes location")
        filetypes = self.tr("csv file (*.csv);;All Files (*)")
        fle = QFileDialog.getOpenFileName(self, title, data, filetypes)[0]
        fle = str(fle)
        if not fle:
            return
        fle = os.path.abspath(fle)

        probe_number = self.mdl.getNumberOfMonitoringPoints()

        lst = self.mdl.ImportProbesFromCSV(fle)

        for idx in range(lst):
            new_name = str(probe_number + idx + 1)
            X, Y, Z = self.mdl.getMonitoringPointCoordinates(new_name)
            self.__insertMonitoringPoint(new_name, X, Y, Z)

            if self.case['salome']:
                self.__salomeHandlerAddMonitoringPoint(probe_number + idx, X, Y, Z)

            idx = idx + 1


    def __salomeHandlerAddMonitoringPoint(self, name, X, Y, Z):
        self.case['probes'].addProbe(str(name), [X, Y, Z])


    def __salomeHandlerDeleteMonitoringPoint(self, l2):
        l2.sort()
        r = len(l2)
        for n in range(r):
            name = str(l2[n])
            self.case['probes'].remove(name)
            for i in range(n, r):
                l2[i] = l2[i] - 1


    @pyqtSlot("QModelIndex")
    def slotSelectedActors(self, idx):
        """
        If salome, hightlights monitoring points in the VTK view.
        """
        log.debug("Current selected row -> %s" % (idx.row(),))
        self.case['probes'].unSelectAll()
        self.case['probes'].select(str(idx.row() + 1))
        for index in self.tableViewPoints.selectionModel().selectedRows():
            name = index.row() + 1
            log.debug("select row -> %s" % name)
            self.case['probes'].select(str(name))


    @pyqtSlot(bool)
    def slotProbesDisplay(self, checked):
        """
        @type checked: C{True} or C{False}
        @param checked: if C{True}, shows the QGroupBox mesh probes display parameters
        """
        if checked:
            self.case['probes'].setVisibility(1)
        else:
            self.case['probes'].setVisibility(0)


    @pyqtSlot('QString')
    def slotProbesRadius(self, text):
        """
        @type text: C{QString}
        @param text: radius for display probes
        """
        if self.lineEditProbesRadius.validator().state == QValidator.Acceptable:
            r = from_qvariant(text, float)
            self.case['probes'].setRadius(r)


    def isSteady(self):
        """
        """
        steady = 0

        if self.case.xmlRootNode().tagName != "NEPTUNE_CFD_GUI" :
            from code_saturne.model.TimeStepModel import TimeStepModel
            idtvar = TimeStepModel(self.case).getTimePassing()
            if idtvar in [-1, 2]:
                steady = 1

        return steady


    @pyqtSlot(int)
    def slotchanged(self, index):
        """
        Changed tab
        """
        self.case['current_tab'] = index


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
