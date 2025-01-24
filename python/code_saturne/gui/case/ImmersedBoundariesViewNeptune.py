# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2024 EDF S.A.
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
This module defines the Immersed boundaries view data management.

This module contains the following classes and function:
- ImmersedBoundariesViewNeptune
"""

#-------------------------------------------------------------------------------
# Standard modules
#-------------------------------------------------------------------------------

import logging, os

#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------

from code_saturne.gui.base.QtCore    import *
from code_saturne.gui.base.QtGui     import *
from code_saturne.gui.base.QtWidgets import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.model.Common import LABEL_LENGTH_MAX, GuiParam
from code_saturne.gui.base.QtPage import IntValidator, DoubleValidator, RegExpValidator, ComboModel
from code_saturne.gui.base.QtPage import from_qvariant, to_text_string
from code_saturne.gui.case.ImmersedBoundariesNeptune import Ui_ImmersedBoundariesNeptune
from code_saturne.model.ImmersedBoundariesModel import ImmersedBoundariesModel
from code_saturne.gui.case.QMegEditorView import QMegEditorView

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("ImmersedBoundariesViewNeptune")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# QLineEdit delegate to attach a label to the solid object
#-------------------------------------------------------------------------------

class ObjectNameDelegate(QItemDelegate):

    def __init__(self, parent = None):
        super(ObjectNameDelegate, self).__init__(parent)
        self.parent = parent


    def createEditor(self, parent, option, index):
        editor = QLineEdit(parent)
        return editor


    def setEditorData(self, editor, index):
        self.value = from_qvariant(index.model().data(index, Qt.DisplayRole),
                                   to_text_string)
        editor.setText(self.value)


    def setModelData(self, editor, model, index):
        value = editor.text()

        if str(value) != "":
            model.setData(index, value, Qt.DisplayRole)


#-------------------------------------------------------------------------------
# QComboBox delegate for the method : Set explicit function or use MED file
#-------------------------------------------------------------------------------

class MethodDelegate(QItemDelegate):
    """
    Use of a combobox to set the method to track the solid
    """

    def __init__(self, parent, mdl):
        super(MethodDelegate, self).__init__(parent)
        self.parent  = parent
        self.mdl     = mdl


    def createEditor(self, parent, option, index):
        editor = QComboBox(parent)

        for itm in ["defined by function",
                    "MEDCoupling using a MED file",
                    "STL file"]:
            editor.addItem(itm)

        editor.installEventFilter(self)
        return editor


    def setEditorData(self, comboBox, index):
        row = index.row()
        col = index.column()
        string = index.model().dataObject[row][col]
        comboBox.setEditText(string)


    def setModelData(self, comboBox, model, index):
        value = comboBox.currentText()
        model.setData(index, value, Qt.DisplayRole)


#-------------------------------------------------------------------------------
# StandarItemModel class
#-------------------------------------------------------------------------------

class StandardItemModel(QStandardItemModel):

    def __init__(self, model, case, tree):
        """
        """
        QStandardItemModel.__init__(self)

        self.headers = [self.tr("Object name"),
                        self.tr("Solid tracking method")]

        self.tooltip = [self.tr("Name of solid object"),
                        self.tr("Method to track the solid")]

        self.setColumnCount(len(self.headers))
        self.dataObject = []
        self.__model = model
        self.case = case
        self.browser = tree


    def data(self, index, role):
        if not index.isValid():
            return None

        row = index.row()
        col = index.column()

        taille = len(self.dataObject)-1
        if row > taille:
            return None

        # Tooltips
        if role == Qt.ToolTipRole:
            return self.tooltip[index.column()]

        elif role == Qt.CheckStateRole:

            if (index.column() == 2):
                data_ibm = self.dataObject[index.row()][index.column()]

                if data_ibm == 'on':
                    return Qt.Checked
                else:
                    return Qt.Unchecked

        # Display
        elif role == Qt.DisplayRole:
            return self.dataObject[index.row()][index.column()]

        elif role == Qt.TextAlignmentRole:
            return Qt.AlignCenter

        return None


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled

        row = index.row()
        col = index.column()

        taille = len(self.dataObject)-1
        if row > taille:
            return Qt.ItemIsEnabled

        if (col == 2): #solve fsi:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsUserCheckable
        else:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable


    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return self.headers[section]
        return None


    def setData(self, index, value, role):
        if not index.isValid():
            return Qt.ItemIsEnabled

        row = index.row()
        col = index.column()

        num = row + 1

        if col < 2:
            self.dataObject[row][col] = str(from_qvariant(value, to_text_string))
            self.__model.setObjectName(num, self.dataObject[row][0])
            self.__model.setObjectMethod(num, self.dataObject[row][1])
            self.browser.configureTree(self.case)

        id1 = self.index(0, 0)
        id2 = self.index(self.rowCount(), 0)
        self.dataChanged.emit(id1, id2)
        return True


    def getData(self, index):
        row = index.row()
        return self.dataObject[row]


    def addItem(self, object_name, object_method):
        """
        Add a row in the table.
        """

        self.dataObject.append([object_name, object_method])
        row = self.rowCount()
        self.setRowCount(row+1)


    def deleteRow(self, row):
        """
        Delete the row in the model
        """
        del self.dataObject[row]
        row = self.rowCount()
        self.setRowCount(row-1)


    def getItem(self, row):
        return self.dataObject[row]


#-------------------------------------------------------------------------------
# QStandardItemModel for monitoring points QTableView
#-------------------------------------------------------------------------------

class StandardItemModelSTLPoints(QStandardItemModel):
    def __init__(self, model, objId):
        """
        """
        QStandardItemModel.__init__(self)

        self.setColumnCount(4)
        self.data_points = []
        self.ibm = model
        self.objId = objId


    def data(self, index, role):
        if not index.isValid():
            return None

        if role == Qt.DisplayRole:

            row = index.row()
            dico = self.data_points[row]

            if index.column() == 0:
                return dico['n']
            elif index.column() == 1:
                return dico['X']
            elif index.column() == 2:
                return dico['Y']
            elif index.column() == 3:
                return dico['Z']
            else:
                return None

        elif role == Qt.TextAlignmentRole:
            return Qt.AlignCenter

        return None


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
                return self.tr("n")
            elif section == 1:
                return self.tr("X")
            elif section == 2:
                return self.tr("Y")
            elif section == 3:
                return self.tr("Z")
        return None


    def setData(self, index, value, role=None):
        row = index.row()
        num_pt = row + 1

        if index.column() == 0:
            n = from_qvariant(value, int)
            self.data_points[row]['n'] = n
        elif index.column() == 1:
            X = from_qvariant(value, float)
            self.data_points[row]['X'] = X
            self.ibm.setSTLObjectSeedPoint(self.objId, num_pt, x=X)
        elif index.column() == 2:
            Y = from_qvariant(value, float)
            self.data_points[row]['Y'] = Y
            self.ibm.setSTLObjectSeedPoint(self.objId, num_pt, y=Y)
        elif index.column() == 3:
            Z = from_qvariant(value, float)
            self.data_points[row]['Z'] = Z
            self.ibm.setSTLObjectSeedPoint(self.objId, num_pt, z=Z)

        #self.dataChanged.emit(index, index)
        id1 = self.index(0, 0)
        id2 = self.index(self.rowCount(), 0)
        self.dataChanged.emit(id1, id2)
        return True


    def insertData(self, num, X, Y, Z):
        """
        Add a new 'item' into the table.
        """
        dico = {}
        dico['n']    = num
        dico['X']    = X
        dico['Y']    = Y
        dico['Z']    = Z

        self.data_points.append(dico)

        row = self.rowCount()
        self.setRowCount(row + 1)


    def replaceData(self, row, label, X, Y, Z):
        """
        Replace value in an existing 'item' into the table.
        """
        self.data_points[row]['n']    = label
        self.data_points[row]['X']    = X
        self.data_points[row]['Y']    = Y
        self.data_points[row]['Z']    = Z

        row = self.rowCount()
        self.setRowCount(row)


    def getLabel(self, index):
        return self.getData(index)['n']


    def getData(self, index):
        row = index.row()
        dico = self.data_points[row]
        return dico


    def deleteData(self, row):
        """
        Delete the row in the model
        """
        del self.data_points[row]
        row = self.rowCount()
        self.setRowCount(row-1)


    def deleteAllData(self):
        """
        Destroy the contents of the list.
        """
        self.data_points = []
        self.setRowCount(0)

#-------------------------------------------------------------------------------
# QItemDelegate for STL seed points points QTableView
#-------------------------------------------------------------------------------

class STLPointDelegate(QItemDelegate):
    def __init__(self, parent=None, case=None, model=None):

        super(STLPointDelegate, self).__init__(parent)
        self.table = parent
        self.case  = case


    def createEditor(self, parent, option, index):
        editor = QLineEdit(parent)
        return editor


    def setEditorData(self, editor, index):
        self.value = from_qvariant(index.model().data(index, Qt.DisplayRole),
                                   to_text_string)
        editor.setText(self.value)


    def setModelData(self, editor, model, index):
        value = editor.text()

        if str(value) != "":
            model.setData(index, value, Qt.DisplayRole)
            dico = model.data_points[index.row()]
            x = float(dico['X'])
            y = float(dico['Y'])
            z = float(dico['Z'])
            label = str(dico['n'])


#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class ImmersedBoundariesViewNeptune(QWidget, Ui_ImmersedBoundariesNeptune):
    """
    """
    def __init__(self, parent, case, stbar, tree):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_ImmersedBoundariesNeptune.__init__(self)
        self.setupUi(self)

        self.case = case
        self.case.undoStopGlobal()
        self.stbar = stbar
        self.ibm = ImmersedBoundariesModel(self.case)
        self.current_obj = None
        self.browser = tree

        # Connections
        self.pushButtonAdd.clicked.connect(self.slotAdd)
        self.pushButtonDelete.clicked.connect(self.slotDelete)
        self.pushButtonExplicit.clicked.connect(self.slotExplicitFormula)
        self.toolButtonMEDFile.clicked.connect(self.slotSearchMEDMesh)
        self.toolButtonSTLFile.clicked.connect(self.slotSearchSTLMesh)
        self.pushButtonAddSTLPoint.clicked.connect(self.slotAddSTLPoint)
        self.pushButtonDeleteSTLPoint.clicked.connect(self.slotDeleteSTLPoint)

        self.dim = ComboModel(self.comboBoxDim, 1, 1)
        self.dim.addItem("3D computation")
        self.dim.addItem("2D computation with symmetry in X-direction")
        self.dim.addItem("2D computation with symmetry in Y-direction")
        self.dim.addItem("2D computation with symmetry in Z-direction")
        self.comboBoxDim.activated[str].connect(self.slotIBMDim)

        self.lineEditSTLFile.setEnabled(False)
        self.lineEditMEDFile.setEnabled(False)

        # Models
        self.model = StandardItemModel(self.ibm, self.case, self.browser)
        self.tableView.setModel(self.model)

        for obj in range(1,self.ibm.getNumberOfObjects()+1):
            self.model.addItem(self.ibm.getObjectName(obj),
                               self.ibm.getObjectMethod(obj))

        if QT_API == "PYQT4":
            self.tableView.verticalHeader().setResizeMode(QHeaderView.ResizeToContents)
            self.tableView.horizontalHeader().setResizeMode(QHeaderView.ResizeToContents)
        elif QT_API == "PYQT5":
            self.tableView.verticalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
            self.tableView.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)

        self.model.dataChanged.connect(self.dataChanged)

        delegateObjectLabel = ObjectNameDelegate(self.tableView)
        self.tableView.setItemDelegateForColumn(0, delegateObjectLabel)

        delegateObjectMethod = MethodDelegate(self.tableView, self.ibm)
        self.tableView.setItemDelegateForColumn(1, delegateObjectMethod)

        self.checkBoxActivate.stateChanged.connect(self.slotCheckActivate)

        self.tableView.clicked[QModelIndex].connect(self.slotChangedSelection)

        # Check for MEDCoupling presence
        self.has_medcoupling = False
        try:
            cfg = self.case.case['package'].config
            self.has_medcoupling = cfg.libs['medcoupling'].have
        except Exception:  # if case/package not available (should not happen)
            print("Warning: package configuration not available")
            pass

        #STL seed points
        self.tableView.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.tableView.selectionModel().selectionChanged.connect(self.on_selection_changed)

        if QT_API == "PYQT4":
            self.tableViewSTLPoints.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        elif QT_API == "PYQT5":
            self.tableViewSTLPoints.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)

        delegateSTLObjectPoints = STLPointDelegate(self.tableViewSTLPoints)
        self.tableViewSTLPoints.setItemDelegate(delegateSTLObjectPoints)

        self.updatePageView()
        self.case.undoStartGlobal()



    def on_selection_changed(self, selected, deselected):
        indexes = self.tableView.selectionModel().selectedRows()

        if indexes:
            obj = indexes[0].row()+1
            self.STLPointsmodel = StandardItemModelSTLPoints(self.ibm, obj)
            self.tableViewSTLPoints.setModel(self.STLPointsmodel)

            #self.STLPointsmodel.deleteAllData()
            for pts_id in range(1,self.ibm.getNumberOfSTLObjectSeedPoints(obj)+1):
                x, y, z = self.ibm.getSTLObjectSeedPoint(obj, pts_id)
                self.STLPointsmodel.insertData(pts_id-1, x, y, z)


    def selectIBMMeshFiles(self, extension):
        """
        Open a File Dialog in order to select mesh files.
        """

        # Meshes directory

        self.mesh_dirs = [None]

        case_dir = self.case['case_path']
        study_path = os.path.split(case_dir)[0]
        d = os.path.join(study_path, 'MESH')
        mesh_files = []

        title = self.tr("Select input IBM mesh file(s)")

        default = self.mesh_dirs[0]
        if default is None:
            default = os.path.split(self.case['case_path'])[0]

        if hasattr(QFileDialog, 'ReadOnly'):
            options  = QFileDialog.DontUseNativeDialog | QFileDialog.ReadOnly
        else:
            options  = QFileDialog.DontUseNativeDialog

        l_mesh_dirs = []
        for i in range(0, len(self.mesh_dirs)):
            if self.mesh_dirs[i] != None:
                l_mesh_dirs.append(QUrl.fromLocalFile(self.mesh_dirs[i]))

        dialog = QFileDialog()
        dialog.setWindowTitle(title)
        dialog.setDirectory(default)

        if extension == 'stl':
            dialog.setNameFilter(self.tr("*.stl"))
        if extension == 'med':
            dialog.setNameFilter(self.tr("*.med"))

        if hasattr(dialog, 'setOptions'):
            dialog.setOptions(options)
        dialog.setSidebarUrls(l_mesh_dirs)
        dialog.setFileMode(QFileDialog.ExistingFiles)

        if dialog.exec_() == 1:
            s = dialog.selectedFiles()
            count = len(s)
            for i in range(count):
                el = str(s[0])
                s = s[1:]
                mesh_files.append(el)

        return mesh_files


    @pyqtSlot("QModelIndex")
    def slotChangedSelection(self, index):
        """
        detect change in selection and update view
        """
        row = self.tableView.currentIndex().row()
        self.current_obj = row + 1

        self.updatePageView()


    @pyqtSlot(int)
    def slotCheckActivate(self, val):

        # Set the method state
        if val == 0:
            self.ibm.setOnOff('off')
        else:
            self.ibm.setOnOff('on')

        # Update the view if needed
        self.updatePageView()


    def dataChanged(self, topLeft, bottomRight):
        self.updatePageView()


    def updatePageView(self):

        if self.ibm.getOnOff() == 'off':
            self.checkBoxActivate.setChecked(False)
            self.groupBoxObjects.hide()
            self.groupBoxExplicit.hide()
            self.groupBoxMEDCoupling.hide()
            self.groupBoxSTL.hide()
            self.comboBoxDim.hide()

        else:
            self.comboBoxDim.show()
            self.checkBoxActivate.setChecked(True)
            self.groupBoxObjects.show()

            #Dimension
            dim = str(self.ibm.getIBMDim())
            self.comboBoxDim.setCurrentText(dim)

            if (self.current_obj is not None
                and self.ibm.getNumberOfObjects() > 0):

                if self.ibm.getObjectMethod(self.current_obj) == 'defined by function':
                    self.groupBoxExplicit.show()
                    self.groupBoxMEDCoupling.hide()
                    self.groupBoxSTL.hide()

                    exp = self.ibm.getObjectFormula(self.current_obj-1)
                    if exp:
                        self.pushButtonExplicit.setToolTip(exp)
                        self.pushButtonExplicit.setStyleSheet("background-color: green")
                    else:
                        self.pushButtonExplicit.setStyleSheet("background-color: red")

                elif self.ibm.getObjectMethod(self.current_obj) == 'MEDCoupling using a MED file':
                    self.groupBoxExplicit.hide()
                    self.groupBoxMEDCoupling.show()
                    self.groupBoxSTL.hide()

                    if self.ibm.getMesh(self.current_obj).endswith('.med'):
                        self.lineEditMEDFile.setText(os.path.basename(
                            self.ibm.getMesh(self.current_obj)))
                    else:
                        self.lineEditMEDFile.setText('')

                elif self.ibm.getObjectMethod(self.current_obj) == 'STL file':
                    self.groupBoxExplicit.hide()
                    self.groupBoxMEDCoupling.hide()
                    self.groupBoxSTL.show()

                    if self.ibm.getMesh(self.current_obj).endswith('.stl'):
                        self.lineEditSTLFile.setText(os.path.basename(
                            self.ibm.getMesh(self.current_obj)))
                    else:
                        self.lineEditSTLFile.setText('')


            else:
                self.groupBoxExplicit.hide()
                self.groupBoxMEDCoupling.hide()
                self.groupBoxSTL.hide()


    @pyqtSlot()
    def slotAdd(self):

        name = self.ibm.defaultValues()['object_name'] \
            + str(self.ibm.getNumberOfObjects()+1)

        #if the name is already taken, change the object_name
        object_list = self.ibm.getObjectsNameList();
        i = 1
        while (name in object_list):
            #name  = '_'.join([self.ibm.defaultValues()['object_name'], str(i)])
            name  = self.ibm.defaultValues()['object_name']+str(i)
            i = i+1

        method = self.ibm.defaultValues()['object_method']

        is_fsi = self.ibm.defaultValues()['object_is_fsi']

        moving_type = self.ibm.defaultValues()['object_moving']

        num = self.ibm.addObject(name, method, is_fsi, moving_type)

        self.model.addItem(name, method)

        self.browser.configureTree(self.case)


    @pyqtSlot()
    def slotDelete(self):
        row = self.tableView.currentIndex().row()

        log.debug("slotDelete -> %s" % (row,))
        if row == -1:
            title = self.tr("Warning")
            msg   = self.tr("You must select an existing object")
            QMessageBox.information(self, title, msg)
        else:
            self.model.deleteRow(row)
            self.ibm.deleteObject(row+1)

        self.browser.configureTree(self.case)


    @pyqtSlot()
    def slotExplicitFormula(self):
        """
        Explicit formula for variable porosity
        """

        objId = self.current_obj

        exp, req, sym = self.ibm.getIBMFormulaComponents(objId-1)
        #exa = """if (x < 0.5)
        #           indicator = 0;
        #         else
        #           indicator = 1;"""

        exa = (
            "# FSI Cylinder with radius r=0.2m in the xy-plane.\n"
            "# indicator values have to be 0 or 1: 1 means only solid and 0 only fluid.\n\n"
            "# To use the FSI center of gravity 'cog_fsi', the FSI option "
            "must be activated in immersed volume conditions page.\n"
            "# Then, 'cog_fsi' must be filled in the FSI tab "
            "associated to the current object.\n\n"
            "r = 0.2;\n"
            "dx = x - cog_fsi_x;\n"
            "dy = y - cog_fsi_y;\n\n"
            "if (dx*dx + dy*dy < r*r)\n"
            "    indicator = 1;"
        )

        name = self.ibm.getObjectName(objId)

        dialog = QMegEditorView(parent        = self,
                                function_type = 'ibm',
                                zone_name     = name,
                                variable_name = 'porosity',
                                expression    = exp,
                                required      = req,
                                symbols       = sym,
                                known_fields  = [],
                                examples      = exa)

        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotExplicitFormula -> %s" % str(result))
            self.pushButtonExplicit.setStyleSheet("background-color: green")
            self.pushButtonExplicit.setToolTip(exp)
            self.ibm.setObjectFormula(objId-1, result)


    @pyqtSlot()
    def slotSearchMEDMesh(self):
        msg = self.tr("Select a mesh file.")
        self.stbar.showMessage(msg, 2000)
        num = self.current_obj

        file_name = self.selectIBMMeshFiles('med')

        if (file_name != [] and file_name[0].endswith('.med')):
            self.ibm.setMesh(num, file_name[0])
            self.lineEditMEDFile.setText(os.path.basename(file_name[0]))


    @pyqtSlot()
    def slotSearchSTLMesh(self):
        msg = self.tr("Select a mesh file.")
        self.stbar.showMessage(msg, 2000)
        num = self.current_obj

        file_name = self.selectIBMMeshFiles('stl')

        if (file_name != []):
            if (file_name[0].endswith('.stl')):
                self.ibm.setMesh(num, file_name[0])
                self.lineEditSTLFile.setText(os.path.basename(file_name[0]))


    @pyqtSlot(str)
    def slotIBMDim(self):
        """
        """
        val = str(self.comboBoxDim.currentText())
        self.ibm.setIBMDim(val)


    @pyqtSlot()
    def slotAddSTLPoint(self):
        """
        Add one STL seed point with these coordinates in the list
        The number of the monitoring point is added at the precedent one
        """

        objId = self.current_obj
        num = self.ibm.getNumberOfSTLObjectSeedPoints(objId)

        self.STLPointsmodel.insertData(num, str('0'), str('0'), str('0'))
        self.ibm.addSTLSeedPoint(objId, num, x=0.0, y=0.0, z=0.0)


    @pyqtSlot()
    def slotDeleteSTLPoint(self):
        """
        Just delete the current selected seed point
        """

        row_pt = self.tableViewSTLPoints.currentIndex().row()
        row_obj = self.tableView.currentIndex().row()

        log.debug("slotDeleteSTLPoint -> %s" % (row_pt,))
        if row_pt == -1:
            title = self.tr("Warning")
            msg   = self.tr("You must select an existing object")
            QMessageBox.information(self, title, msg)
        else:
            self.STLPointsmodel.deleteData(row_pt)
            self.ibm.deleteSTLSeedPoint(row_obj+1,row_pt+1)


#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
