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
This module contains the following classes:
- ProfilesView
"""

#-------------------------------------------------------------------------------
# Standard modules
#-------------------------------------------------------------------------------

import os, sys
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

from code_saturne.model.Common import LABEL_LENGTH_MAX, GuiParam
from code_saturne.gui.base.QtPage import IntValidator, DoubleValidator, RegExpValidator, ComboModel
from code_saturne.gui.case.ProfilesForm import Ui_ProfilesForm
from code_saturne.model.ProfilesModel import ProfilesModel
from code_saturne.gui.case.QMegEditorView import QMegEditorView
from code_saturne.model.NotebookModel import NotebookModel
from code_saturne.model.PyVistaWrapper import pyvistaWrapper

if QT_API == "PYQT6":
    from code_saturne.gui.case import resources_pages_rc

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("ProfilesView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# StandarItemModel class
#-------------------------------------------------------------------------------

class StandardItemModelProfile(QStandardItemModel):

    def __init__(self):
        """
        """
        QStandardItemModel.__init__(self)
        self.setColumnCount(2)
        self.dataProfile = []


    def data(self, index, role):
        if not index.isValid():
            return None
        if role == Qt.ItemDataRole.DisplayRole:
            return self.dataProfile[index.row()][index.column()]
        elif role == Qt.ItemDataRole.TextAlignmentRole:
            return Qt.AlignmentFlag.AlignCenter
        return None


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemFlag.ItemIsEnabled
        else:
            return Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsSelectable


    def headerData(self, section, orientation, role):
        if orientation == Qt.Orientation.Horizontal and role == Qt.ItemDataRole.DisplayRole:
            if section == 0:
                return self.tr("Filename")
            elif section == 1:
                return self.tr("Variables")
        elif role == Qt.ItemDataRole.TextAlignmentRole:
            return Qt.AlignmentFlag.AlignCenter
        return None


    def setData(self, index, value, role):
        self.dataChanged.emit(index, index)
        return True


    def addItem(self, label, prof):
        """
        Add a row in the table.
        """
        self.dataProfile.append([label, prof])
        row = self.rowCount()
        self.setRowCount(row+1)


    def replaceItem(self, row, label, prof):
        """
        Replace a row in the table.
        """
        self.dataProfile[row] = [label, prof]


    def deleteRow(self, row):
        """
        Delete the row in the model
        """
        del self.dataProfile[row]
        row = self.rowCount()
        self.setRowCount(row-1)


    def getItem(self, row):
        """
        Returns the name of the mesh file.
        """
        return self.dataProfile[row]


    def getLabel(self, row):
        """
        Returns the name of the mesh file.
        """
        return self.dataProfile[row][0]

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class ProfilesView(QWidget, Ui_ProfilesForm):
    """
    """
    def __init__(self, parent, case, stbar):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_ProfilesForm.__init__(self)
        self.setupUi(self)

        self.case = case
        self.case.undoStopGlobal()
        self.mdl = ProfilesModel(self.case)
        self.notebook = NotebookModel(self.case)

        #  Initialize variables concerning the display of the Hlist

        self.entriesNumber = 0

        # Models
        self.modelProfile = StandardItemModelProfile()
        self.treeViewProfile.setModel(self.modelProfile)
        self.treeViewProfile.resizeColumnToContents(0)

        # QListView layout
        self.gridlayout1 = QGridLayout(self.widgetDrag)
        self.gridlayout1.setContentsMargins(0, 0, 0, 0)
        self.DragList = QListView(self.widgetDrag)
        self.gridlayout1.addWidget(self.DragList,0,0,1,1)

        self.gridlayout2 = QGridLayout(self.widgetDrop)
        self.gridlayout2.setContentsMargins(0, 0, 0, 0)
        self.DropList = QListView(self.widgetDrop)
        self.gridlayout2.addWidget(self.DropList,0,0,1,1)

        self.modelDrag = QStringListModel()
        self.modelDrop = QStringListModel()
        self.DragList.setModel(self.modelDrag)
        self.DropList.setModel(self.modelDrop)
        self.DragList.setAlternatingRowColors(True)
        self.DragList.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.DropList.setAlternatingRowColors(True)
        self.DropList.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)

        # Combo items
        self.modelFreq = ComboModel(self.comboBoxFreq, 3, 1)
        self.modelFreq.addItem(self.tr("at the end of the calculation"), "end")
        self.modelFreq.addItem(self.tr("at each 'n' time steps"), "frequency")
        self.modelFreq.addItem(self.tr("Output every 'x' seconds"), 'time_value')

        self.modelFormat = ComboModel(self.comboBoxFormat, 2, 1)
        self.modelFormat.addItem(self.tr(".dat"), "DAT")
        self.modelFormat.addItem(self.tr(".csv"), "CSV")

        # Connections
        self.treeViewProfile.pressed[QModelIndex].connect(self.slotSelectProfile)
        self.pushButtonAdd.clicked.connect(self.slotAddProfile)
        self.pushButtonDelete.clicked.connect(self.slotDeleteProfile)
        self.pushButtonAddVar.clicked.connect(self.slotAddVarProfile)
        self.pushButtonSuppressVar.clicked.connect(self.slotDeleteVarProfile)
        self.comboBoxFreq.activated[int].connect(self.slotFrequencyType)
        self.comboBoxFormat.activated[int].connect(self.slotFormatType)
        self.pushButtonFormula.clicked.connect(self.slotFormula)

        # Mesh viewer button and display mode
        self.pyvista = pyvistaWrapper()
        self.pushButtonVisuProfile = QPushButton(self.tr("View profile"))
        self.pushButtonVisuProfile.setToolTip(
            self.tr("View the 1D profile inside the 3D mesh volume."))
        self.pushButtonVisuProfile.setEnabled(self.pyvista.isActive())
        self.pushButtonVisuProfile.clicked.connect(self.slotVisuProfile)
        self.horizontalLayout.addWidget(self.pushButtonVisuProfile)

        self.comboBoxProfileDisplay = QComboBox()
        self.comboBoxProfileDisplay.addItem(self.tr("Straight line"), "line")
        self.comboBoxProfileDisplay.addItem(self.tr("Points"), "points")
        self.comboBoxProfileDisplay.setToolTip(
            self.tr("Choose profile representation: continuous straight line or discrete probe points."))
        self.comboBoxProfileDisplay.setEnabled(self.pyvista.isActive())
        self.horizontalLayout.addWidget(self.comboBoxProfileDisplay)

        self.lineEditBaseName.textChanged[str].connect(self.slotBaseName)
        self.lineEditFreq.textChanged[str].connect(self.slotFrequence)
        self.lineEditFreqTime.textChanged[str].connect(self.slotFrequenceTime)
        self.lineEditNbPoint.textChanged[str].connect(self.slotNbPoint)

        self.checkBoxSnapToCenter.clicked.connect(self.slotSnapToCenter)
        self.checkBoxActivateInterpolation.clicked.connect(self.slotActivateInterpolation)

        # Validators
        validatorFreq = IntValidator(self.lineEditFreq, min=0)
        validatorFreq.setExclusiveMin(True)
        self.lineEditFreq.setValidator(validatorFreq)

        validatorFreqT = DoubleValidator(self.lineEditFreqTime, min=0.)
        validatorFreqT.setExclusiveMin(True)
        self.lineEditFreqTime.setValidator(validatorFreqT)

        validatorNbPoint = IntValidator(self.lineEditNbPoint, min=0)
        self.lineEditNbPoint.setValidator(validatorNbPoint)

        rx = r"[\-_A-Za-z0-9]{1," + str(LABEL_LENGTH_MAX) + "}"
        validatorBaseName =  RegExpValidator(self.lineEditBaseName, QRegularExpression(rx))
        self.lineEditBaseName.setValidator(validatorBaseName)

        #update list of variables, properties, scalars ...
        liste_label = []
        for label in self.mdl.getVariablesAndVolumeProperties():
            liste_label.append(label)
        self.modelDrag.setStringList(sorted(liste_label, key=str.lower))

        #update list of profiles for view from xml file
        for lab in self.mdl.getProfilesLabelsList():
            self.entriesNumber = self.entriesNumber + 1
            label, fmt, lst, choice, freq, formula, nb_point = self.mdl.getProfileData(lab)
            self.__insertProfile(label, lst)

        self.__eraseEntries()

        self.case.undoStartGlobal()


    def __verifLabel(self):
        """
        Verif label.
        """
        label = str(self.lineEditBaseName.text())
        if label in self.mdl.getProfilesLabelsList():
            default = {}
            default['label'] = label
            default['list'] = self.mdl.getProfilesLabelsList()
            rx = r"[\-_A-Za-z0-9]{1," + str(LABEL_LENGTH_MAX) + "}"
            default['regexp'] = QRegularExpression(rx)
            from code_saturne.gui.case.VerifyExistenceLabelDialogView import VerifyExistenceLabelDialogView
            dialog = VerifyExistenceLabelDialogView(self, default)
            if dialog.exec():
                result = dialog.get_result()
                label = result['label']
                if result['label'] == default['label']:
                    label = ""
        return label


    @Slot(int)
    def slotFrequencyType(self, idx):
        """
        Input choice for frequency for profile.
        """
        text = self.comboBoxFreq.currentText()
        choice = self.modelFreq.dicoV2M[str(text)]

        if choice == "end":
            nfreq = -1
            self.lineEditFreq.setText(str(nfreq))
            self.lineEditFreq.show()
            self.lineEditFreqTime.hide()
            self.lineEditFreq.setDisabled(True)

        elif choice == "frequency":
            self.lineEditFreq.show()
            self.lineEditFreqTime.hide()
            nfreq = self.mdl.getOutputFrequency(self.label_select)
            if nfreq == -1:
                nfreq = 1
            self.lineEditFreq.setEnabled(True)
            self.lineEditFreq.setText(str(nfreq))

        elif choice == "time_value":
            self.lineEditFreq.hide()
            self.lineEditFreqTime.show()
            nfreq = self.mdl.getOutputFrequency(self.label_select)
            if nfreq == -1:
                nfreq = 1.0
            self.lineEditFreqTime.setText(str(nfreq))

        self.mdl.setOutputType(self.label_select, choice)
        self.mdl.setOutputFrequency(self.label_select, nfreq)


    @Slot(int)
    def slotFormatType(self, idx):
        """
        Input choice for frequency for profile.
        """
        text = self.comboBoxFormat.currentText()
        fmt = self.modelFormat.dicoV2M[str(text)]
        self.mdl.setFormat(self.label_select, fmt)


    def __infoProfile(self, row):
        """
        Return info from the argument entry.
        """
        label = self.modelProfile.getLabel(row)
        lab, fmt, lst, choice, freq, formula, nb_point = self.mdl.getProfileData(label)
        return label, fmt, lst, choice, freq, formula, nb_point


    def __insertProfile(self, label, lst):
        """
        Insert values in table view.
        """
        self.modelProfile.addItem(label, " ; ".join(lst))


    @Slot()
    def slotAddProfile(self):
        """
        Set in view label and variables to see on profile
        """
        var_prof = []
        label = self.mdl.addProfile()
        self.__insertProfile(label, var_prof)
        self.__eraseEntries()


    @Slot()
    def slotDeleteProfile(self):
        """
        Delete the profile from the list (one by one).
        """
        row = self.treeViewProfile.currentIndex().row()
        log.debug("slotDeleteProfile -> %s" % (row,))
        if row == -1:
            title = self.tr("Warning")
            msg   = self.tr("You must select an existing profile")
            QMessageBox.information(self, title, msg)
        else:
            label, fmt, lst, choice, freq, formula, nb_point = self.__infoProfile(row)
            self.modelProfile.deleteRow(row)
            self.mdl.deleteProfile(label)
            self.__eraseEntries()


    @Slot(QModelIndex)
    def slotSelectProfile(self, index):
        """
        Return the selected item from the list.
        """
        self.groupBoxProfile.show()

        row = index.row()
        log.debug("slotSelectProfile -> %s" % (row,))

        label, fmt, liste, choice, freq, formula, nb_point = self.__infoProfile(row)
        self.label_select = label

        self.lineEditBaseName.setText(str(label))
        self.modelFormat.setItem(str_model=fmt)

        self.modelFreq.setItem(str_model=choice)
        if choice == "end":
            self.lineEditFreq.show()
            self.lineEditFreqTime.hide()
            self.lineEditFreq.setText(str("-1"))
            self.lineEditFreq.setDisabled(True)

        elif choice == "frequency":
            self.lineEditFreq.show()
            self.lineEditFreqTime.hide()
            self.lineEditFreq.setEnabled(True)
            self.lineEditFreq.setText(str(freq))

        elif choice == "time_value":
            self.lineEditFreq.hide()
            self.lineEditFreqTime.show()
            self.lineEditFreqTime.setText(str(freq))

        self.lineEditNbPoint.setText(str(nb_point))

        if self.mdl.getSnapMode(label) == "snap_to_center":
            self.checkBoxSnapToCenter.setChecked(True)
        else:
            self.checkBoxSnapToCenter.setChecked(False)

        if self.mdl.getProfileInterpolation(label) == "yes":
            self.checkBoxActivateInterpolation.setChecked(True)
        else:
            self.checkBoxActivateInterpolation.setChecked(False)


        self.modelDrop.setStringList([])
        liste = [str(s) for s in liste]

        self.modelDrop.setStringList(liste)

        exp = self.mdl.getFormula(self.label_select)
        if exp:
            self.pushButtonFormula.setStyleSheet("background-color: green")
            self.pushButtonFormula.setToolTip(exp)
        else:
            self.pushButtonFormula.setStyleSheet("background-color: red")

        self.pushButtonVisuProfile.setEnabled(self.pyvista.isActive())
        if hasattr(self, 'comboBoxProfileDisplay'):
            self.comboBoxProfileDisplay.setEnabled(self.pyvista.isActive())


    @Slot()
    def slotAddVarProfile(self):
        """
        Add a new var from list to profile
        """
        if (self.DragList.currentIndex().row() >=0) :
            liste = self.modelDrop.stringList()
            var = self.modelDrag.stringList()[self.DragList.currentIndex().row()]
            if var not in liste :
                liste.append(var)
            liste = [str(s) for s in liste]
            self.modelDrop.setStringList(liste)
            self.mdl.setVariable(self.label_select, liste)

            row = self.treeViewProfile.currentIndex().row()
            liste = self.mdl.getVariable(self.label_select)
            self.modelProfile.replaceItem(row, self.label_select, " ; ".join(liste))


    @Slot()
    def slotDeleteVarProfile(self):
        """
        Supress a var from profile
        """
        self.modelDrop.removeRows(self.DropList.currentIndex().row(), 1)
        liste = self.modelDrop.stringList()
        liste = [str(s) for s in liste]
        self.mdl.setVariable(self.label_select, liste)

        row = self.treeViewProfile.currentIndex().row()
        liste = self.mdl.getVariable(self.label_select)
        self.modelProfile.replaceItem(row, self.label_select, " ; ".join(liste))


    def __eraseEntries(self):
        """
        Delete all caracters in the entries.
        """
        self.groupBoxProfile.hide()
        self.label_select = None
        self.treeViewProfile.clearSelection()


    @Slot()
    def slotFormula(self):
        """
        """
        exp = self.mdl.getFormula(self.label_select)
        exa = """#example: a line segment between (x0, y0, z0) and (x1, y1, z1)
#(s, the parameter is always between 0 and 1)
x0 = 3.2;
x1 = 2;
y0 = 2;
y1 = 2;
z0 = 5;
z1 = -0.5;
x = x1*s + x0*(1.-s);
y = y1*s + y0*(1.-s);
z = z1*s + z0*(1.-s);"""
        req = [('x', "x formula"),
               ('y', "y formula"),
               ('z', "z formula")]
        sym = [('s', 'parameter')]

        for (nme, val) in self.notebook.getNotebookList():
            sym.append((nme, 'value (notebook) = ' + str(val)))

        dialog = QMegEditorView(self,
                                function_type = 'pfl',
                                zone_name     = self.label_select,
                                variable_name = 'coords',
                                expression = exp,
                                required   = req,
                                symbols    = sym,
                                examples   = exa)
        if dialog.exec():
            result = dialog.get_result()
            log.debug("slotLineFormula -> %s" % str(result))
            self.mdl.setFormula(self.label_select, str(result))
            self.pushButtonFormula.setToolTip(result)
            self.pushButtonFormula.setStyleSheet("background-color: green")


    @Slot(str)
    def slotBaseName(self, text):
        """
        """
        lst = self.mdl.getProfilesLabelsList()
        if text not in lst:
            if self.lineEditBaseName.validator().state == QValidator.State.Acceptable:
                self.mdl.setLabel(self.label_select, str(text))
                self.label_select = str(text)

                row = self.treeViewProfile.currentIndex().row()
                liste = self.mdl.getVariable(self.label_select)
                self.modelProfile.replaceItem(row, self.label_select, " ; ".join(liste))


    @Slot(str)
    def slotFrequence(self, text):
        """
        """
        if self.lineEditFreq.validator().state == QValidator.State.Acceptable:
            self.mdl.setOutputFrequency(self.label_select, int(text))


    @Slot(str)
    def slotFrequenceTime(self, text):
        """
        """
        if self.lineEditFreqTime.validator().state == QValidator.State.Acceptable:
            self.mdl.setOutputFrequency(self.label_select, float(text))


    @Slot(str)
    def slotNbPoint(self, text):
        """
        """
        if self.lineEditNbPoint.validator().state == QValidator.State.Acceptable:
            self.mdl.setNbPoint(self.label_select, int(text))


    @Slot()
    def slotSnapToCenter(self):
        """
        """
        if self.checkBoxSnapToCenter.isChecked():
            self.mdl.setSnapMode(self.label_select, "snap_to_center")
        else:
            self.mdl.setSnapMode(self.label_select, "none")

    @Slot()
    def slotActivateInterpolation(self):
        """
        """
        if self.checkBoxActivateInterpolation.isChecked():
            self.mdl.setProfileInterpolation(self.label_select, "yes")
        else:
            self.mdl.setProfileInterpolation(self.label_select, "no")

    @Slot()
    def slotVisuProfile(self, checked=False, display_mode=None):
        """
        Open the 3D mesh viewer to display the defined 1D profile inside the mesh volume.
        """
        import os
        from code_saturne.gui.case.MeshViewerDialog import MeshViewerDialog, eval_profile_formula
        from code_saturne.model.SolutionDomainModel import SolutionDomainModel

        if not isinstance(display_mode, str):
            display_mode = None
        if display_mode is None and hasattr(self, 'comboBoxProfileDisplay'):
            display_mode = self.comboBoxProfileDisplay.currentData()
        if not display_mode:
            display_mode = "line"

        if not hasattr(self, 'label_select') or not self.label_select:
            row = self.treeViewProfile.currentIndex().row()
            if row >= 0:
                self.label_select = self.modelProfile.getLabel(row)
            else:
                QMessageBox.information(self, self.tr("Visualization"),
                    self.tr("Please first select a 1D profile from the list."))
                return

        # Retrieve profile data
        lab, fmt, lst, choice, freq, formula, nb_point = self.mdl.getProfileData(self.label_select)

        if not formula or not str(formula).strip():
            QMessageBox.information(self, self.tr("Visualization"),
                self.tr("No mathematical formula defined for profile '%s'.\n"
                        "Define the formula (x, y, z as a function of s) using "
                        "the mathematical expression editor first.") % self.label_select)
            return

        # Notebook variables
        notebook_vars = dict(self.notebook.getNotebookList())

        # Evaluate the active profile
        try:
            points = eval_profile_formula(formula, nb_points=nb_point, notebook_vars=notebook_vars)
        except Exception as e:
            QMessageBox.warning(self, self.tr("Formula Evaluation Error"),
                self.tr("Error evaluating mathematical formula for profile '%s':\n%s") % (self.label_select, str(e)))
            return

        active_prof = {
            'label': self.label_select,
            'formula': formula,
            'nb_points': nb_point,
            'points': points
        }

        # Also collect other defined profiles if available
        profiles_list = [active_prof]
        for lab in self.mdl.getProfilesLabelsList():
            if lab == self.label_select:
                continue
            l, f, ls, c, fr, form, nbp = self.mdl.getProfileData(lab)
            if form and str(form).strip():
                try:
                    pts = eval_profile_formula(form, nb_points=nbp, notebook_vars=notebook_vars)
                    profiles_list.append({
                        'label': lab,
                        'formula': form,
                        'nb_points': nbp,
                        'points': pts
                    })
                except Exception:
                    pass

        # Reconstruct .csm file path
        mdl = SolutionDomainModel(self.case)
        mesh_list = mdl.getMeshList()

        if not mesh_list:
            QMessageBox.warning(self, self.tr("Mesh not found"),
                self.tr("Case has no associated mesh."))
            return

        mesh_name, mesh_dir = mesh_list[0]
        if mesh_dir:
            mesh_file = os.path.join(mesh_dir, mesh_name)
        else:
            case_path = self.case['case_path']
            mesh_file = os.path.join(case_path, '..', 'MESH', mesh_name)
            mesh_file = os.path.normpath(mesh_file)
            if not os.path.isfile(mesh_file):
                mesh_file = os.path.join(case_path, mesh_name)

        if not os.path.isfile(mesh_file):
            QMessageBox.warning(self, self.tr("File not found"),
                self.tr("Impossible to find :\n%s" % mesh_file))
            return

        if not mesh_file.endswith(".csm"):
            QMessageBox.information(self, self.tr("Unsupported file format"),
                self.tr("Mesh view functionality is only available "
                        "for .csm files."))
            return

        # Open in profile mode
        self._visu_dlg = MeshViewerDialog(self, mesh_file,
                                           profiles=profiles_list,
                                           case=self.case,
                                           profile_display_mode=display_mode)
        self._visu_dlg.exec_()


#-------------------------------------------------------------------------------
# Testing part
#-------------------------------------------------------------------------------

import unittest


class ProfilesViewTestCase(unittest.TestCase):
    """
    Unit test cases for ProfilesView and 1D profile visualization.
    """
    @classmethod
    def setUpClass(cls):
        os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
        cls.app = QApplication.instance() or QApplication([])

    def setUp(self):
        from code_saturne.model.XMLengine import Case
        from code_saturne.model.XMLinitialize import XMLinit
        self.case = Case(None)
        XMLinit(self.case).initialize()
        self.view = ProfilesView(None, self.case, None)

    def tearDown(self):
        if hasattr(self, 'view') and self.view:
            self.view.close()
        del self.view
        del self.case

    def testVisuProfileButtonInitialization(self):
        """
        Verify that pushButtonVisuProfile is properly added and configured.
        """
        self.assertTrue(hasattr(self.view, 'pushButtonVisuProfile'))
        self.assertEqual(self.view.pushButtonVisuProfile.text(), "View profile")
        self.assertTrue(bool(self.view.pushButtonVisuProfile.toolTip()))
        from code_saturne.model.PyVistaWrapper import pyvistaWrapper
        pv = pyvistaWrapper()
        self.assertEqual(self.view.pushButtonVisuProfile.isEnabled(), pv.isActive())

        # Test comboBoxProfileDisplay
        self.assertTrue(hasattr(self.view, 'comboBoxProfileDisplay'))
        self.assertEqual(self.view.comboBoxProfileDisplay.count(), 2)
        self.assertEqual(self.view.comboBoxProfileDisplay.itemText(0), "Straight line")
        self.assertEqual(self.view.comboBoxProfileDisplay.itemData(0), "line")
        self.assertEqual(self.view.comboBoxProfileDisplay.itemText(1), "Points")
        self.assertEqual(self.view.comboBoxProfileDisplay.itemData(1), "points")
        self.assertEqual(self.view.comboBoxProfileDisplay.isEnabled(), pv.isActive())

    def testSlotVisuProfileNoSelection(self):
        """
        Verify QMessageBox when slotVisuProfile is called without a selected profile.
        """
        info_boxes = []
        orig_info = QMessageBox.information
        QMessageBox.information = lambda parent, title, msg: info_boxes.append((title, msg))
        try:
            self.view.slotVisuProfile()
            self.assertEqual(len(info_boxes), 1)
            self.assertIn("select", info_boxes[0][1].lower())
        finally:
            QMessageBox.information = orig_info

    def testSlotVisuProfileNoFormula(self):
        """
        Verify QMessageBox when profile has no mathematical formula defined.
        """
        label = self.view.mdl.addProfile()
        self.view.label_select = label
        orig_get_data = self.view.mdl.getProfileData
        self.view.mdl.getProfileData = lambda lab: (lab, 'CSV', [], 'frequency', 1, '', 10)

        info_boxes = []
        orig_info = QMessageBox.information
        QMessageBox.information = lambda parent, title, msg: info_boxes.append((title, msg))
        try:
            self.view.slotVisuProfile()
            self.assertEqual(len(info_boxes), 1)
            self.assertIn("mathematical formula", info_boxes[0][1].lower())
        finally:
            self.view.mdl.getProfileData = orig_get_data
            QMessageBox.information = orig_info

    def testSlotVisuProfileInvalidFormula(self):
        """
        Verify QMessageBox when profile formula has invalid syntax or variables.
        """
        label = self.view.mdl.addProfile()
        self.view.mdl.setFormula(label, "x = invalid_token * s; y = s; z = s;")
        self.view.label_select = label

        warning_boxes = []
        orig_warning = QMessageBox.warning
        QMessageBox.warning = lambda parent, title, msg: warning_boxes.append((title, msg))
        try:
            self.view.slotVisuProfile()
            self.assertEqual(len(warning_boxes), 1)
            self.assertIn("error evaluating", warning_boxes[0][1].lower())
        finally:
            QMessageBox.warning = orig_warning

    def testSlotVisuProfileNoMesh(self):
        """
        Verify QMessageBox when case has no associated mesh.
        """
        label = self.view.mdl.addProfile()
        self.view.mdl.setFormula(label, "x = s; y = 2*s; z = 3*s;")
        self.view.label_select = label

        warning_boxes = []
        orig_warning = QMessageBox.warning
        QMessageBox.warning = lambda parent, title, msg: warning_boxes.append((title, msg))
        try:
            self.view.slotVisuProfile()
            self.assertEqual(len(warning_boxes), 1)
            self.assertIn("mesh", warning_boxes[0][1].lower())
        finally:
            QMessageBox.warning = orig_warning

    def testSlotVisuProfileExecutionWithMesh(self):
        """
        Verify slotVisuProfile successfully creates and configures MeshViewerDialog.
        """
        csm_file = '/home/i76777/Codes/TEST/2026/PYVISTA/CASE1/RESU/20260520-1104/checkpoint/mesh_input.csm'
        if not os.path.isfile(csm_file):
            self.skipTest("Sample CSM mesh not found at %s" % csm_file)

        from code_saturne.model.SolutionDomainModel import SolutionDomainModel
        from code_saturne.gui.case.MeshViewerDialog import MeshViewerDialog

        # Set mesh on case
        sdm = SolutionDomainModel(self.case)
        sdm.addMesh((os.path.basename(csm_file), os.path.dirname(csm_file)))

        # Add profile
        label = self.view.mdl.addProfile()
        formula = "x = 0.1 * s; y = 0.2 * s; z = 0.3 * s;"
        self.view.mdl.setFormula(label, formula)
        self.view.mdl.setNbPoint(label, 15)
        self.view.label_select = label

        # Also add a second profile to verify multiple profile handling
        label2 = self.view.mdl.addProfile()
        self.view.mdl.setFormula(label2, "x = 0.5 * s; y = 0.5; z = 0.1 * s;")
        self.view.mdl.setNbPoint(label2, 10)

        # Mock dialog exec to avoid blocking modal loop
        exec_called = []
        orig_exec = MeshViewerDialog.exec_
        def mock_exec(dlg_self):
            exec_called.append(dlg_self)
            return 0
        MeshViewerDialog.exec_ = mock_exec

        try:
            self.view.slotVisuProfile()
        finally:
            MeshViewerDialog.exec_ = orig_exec

        self.assertEqual(len(exec_called), 1)
        dlg = exec_called[0]
        self.assertEqual(len(dlg.profiles), 2)
        self.assertEqual(dlg.profiles[0]['label'], label)
        self.assertEqual(len(dlg.profiles[0]['points']), 15)
        self.assertEqual(dlg.profiles[1]['label'], label2)
        self.assertEqual(len(dlg.profiles[1]['points']), 10)
        self.assertEqual(dlg.profile_display_mode, 'line')
        self.assertTrue(dlg.show_wireframe)
        self.assertTrue(dlg.show_volume)
        dlg.close()

    def testSlotVisuProfilePointsModeExecutionWithMesh(self):
        """
        Verify slotVisuProfile forwards profile_display_mode='points' to MeshViewerDialog.
        """
        csm_file = '/home/i76777/Codes/TEST/2026/PYVISTA/CASE1/RESU/20260520-1104/checkpoint/mesh_input.csm'
        if not os.path.isfile(csm_file):
            self.skipTest("Sample CSM mesh not found at %s" % csm_file)

        from code_saturne.model.SolutionDomainModel import SolutionDomainModel
        from code_saturne.gui.case.MeshViewerDialog import MeshViewerDialog

        # Set mesh on case
        sdm = SolutionDomainModel(self.case)
        sdm.addMesh((os.path.basename(csm_file), os.path.dirname(csm_file)))

        # Add profile with 25 probe points
        label = self.view.mdl.addProfile()
        formula = "x = 0.1 * s; y = 0.2 * s; z = 0.3 * s;"
        self.view.mdl.setFormula(label, formula)
        self.view.mdl.setNbPoint(label, 25)
        self.view.label_select = label

        # Select points mode in comboBoxProfileDisplay
        self.view.comboBoxProfileDisplay.setCurrentIndex(1)
        self.assertEqual(self.view.comboBoxProfileDisplay.currentData(), 'points')

        exec_called = []
        orig_exec = MeshViewerDialog.exec_
        def mock_exec(dlg_self):
            exec_called.append(dlg_self)
            return 0
        MeshViewerDialog.exec_ = mock_exec

        try:
            self.view.slotVisuProfile()
        finally:
            MeshViewerDialog.exec_ = orig_exec

        self.assertEqual(len(exec_called), 1)
        dlg = exec_called[0]
        self.assertEqual(dlg.profile_display_mode, 'points')
        self.assertEqual(len(dlg.profiles[0]['points']), 25)
        self.assertTrue(dlg.show_wireframe)
        self.assertTrue(dlg.show_volume)
        dlg.close()


def suite():
    testSuite = unittest.TestSuite()
    testSuite.addTest(unittest.makeSuite(ProfilesViewTestCase))
    return testSuite


def runTest():
    runner = unittest.TextTestRunner()
    runner.run(suite())


if __name__ == "__main__":
    runTest()


#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
