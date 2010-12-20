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
from OutputControlForm import Ui_OutputControlForm
import Base.QtPage as QtPage
from Pages.OutputControlModel import OutputControlModel
from Pages.ConjugateHeatTransferModel import ConjugateHeatTransferModel
from Pages.MobileMeshModel import MobileMeshModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("OutputControlView")
log.setLevel(GuiParam.DEBUG)

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
        self.modelPostProcessing = QtPage.ComboModel(self.comboBoxPostProcessing,2,1)
        self.modelMeshes         = QtPage.ComboModel(self.comboBoxMeshes,3,1)
        self.modelFMTCHR         = QtPage.ComboModel(self.comboBoxFMTCHR,3,1)
        self.modelFormat         = QtPage.ComboModel(self.comboBoxFormat,2,1)
        self.modelPolygon        = QtPage.ComboModel(self.comboBoxPolygon,3,1)
        self.modelPolyhedra      = QtPage.ComboModel(self.comboBoxPolyhedra,3,1)
        self.modelHisto          = QtPage.ComboModel(self.comboBoxHisto,3,1)

        self.modelOutput.addItem(self.tr("No output"), 'None')
        self.modelOutput.addItem(self.tr("Output listing at each time step"), 'At each step')
        self.modelOutput.addItem(self.tr("Output every 'n' time steps"), 'Frequency_l')

        self.modelPostProcessing.addItem(self.tr("Only at the end of calculation"), 'At the end')
        self.modelPostProcessing.addItem(self.tr("At each time step"), 'At each step')
        self.modelPostProcessing.addItem(self.tr("Post-processing every 'n' time steps"), 'Frequency_c')
        self.modelPostProcessing.addItem(self.tr("Post-processing every 'x' second(s)"), 'Frequency_c_x')

        self.modelMeshes.addItem(self.tr("fixed"), '0')
        self.modelMeshes.addItem(self.tr("deformable"), '1')
        self.modelMeshes.addItem(self.tr("modifiable"), '2')
        self.modelMeshes.addItem(self.tr("fixed (with displacement)"), '10')
        self.modelMeshes.addItem(self.tr("deformable (with displacement)"), '11')
        self.modelMeshes.addItem(self.tr("modifiable (with displacement)"), '12')

        ale = self.ale()

        self.modelFMTCHR.addItem(self.tr("EnSight Gold"), 'EnSight')
        self.modelFMTCHR.addItem(self.tr("MED"), 'MED_fichier')
        self.modelFMTCHR.addItem(self.tr("CGNS"), 'CGNS')

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

        # Connections

        self.connect(self.comboBoxOutput, SIGNAL("activated(const QString&)"), self.slotOutputListing)
        self.connect(self.lineEditNTLIST, SIGNAL("textChanged(const QString &)"), self.slotListingFrequency)
        self.connect(self.comboBoxPostProcessing, SIGNAL("activated(const QString&)"), self.slotOutputPostpro)
        self.connect(self.lineEditNTCHR, SIGNAL("textChanged(const QString &)"), self.slotPostproFrequency)
        self.connect(self.lineEditFRCHR, SIGNAL("textChanged(const QString &)"), self.slotPostproFrequencyTime)
        self.connect(self.checkBoxICHRVL, SIGNAL("clicked()"), self.slotVolumeOuput)
        self.connect(self.checkBoxICHRBO, SIGNAL("clicked()"), self.slotBoundaryOuput)
        self.connect(self.checkBoxICHRSY, SIGNAL("clicked()"), self.slotSyrthesOutput)
        self.connect(self.comboBoxMeshes, SIGNAL("activated(const QString&)"), self.slotTypePostMeshes)
        self.connect(self.comboBoxFMTCHR, SIGNAL("activated(const QString&)"), self.slotOutputFormat)
        self.connect(self.comboBoxFormat, SIGNAL("activated(const QString&)"), self.slotOutputOptions)
        self.connect(self.comboBoxPolygon, SIGNAL("activated(const QString&)"), self.slotOutputOptions)
        self.connect(self.comboBoxPolyhedra, SIGNAL("activated(const QString&)"), self.slotOutputOptions)
        self.connect(self.checkBoxBigEndian, SIGNAL("clicked()"), self.slotOutputOptions)


        self.connect(self.toolButtonAdd, SIGNAL("clicked()"), self.slotAddMonitoringPoint)
        self.connect(self.toolButtonDelete, SIGNAL("clicked()"), self.slotDeleteMonitoringPoints)
        self.connect(self.comboBoxHisto, SIGNAL("activated(const QString&)"), self.slotMonitoringPoint)
        self.connect(self.lineEditHisto, SIGNAL("textChanged(const QString &)"), self.slotMonitoringPointFrequency)
        self.connect(self.lineEditFRHisto, SIGNAL("textChanged(const QString &)"), self.slotMonitoringPointFrequencyTime)

        # Validators

        validatorNTLIST = QtPage.IntValidator(self.lineEditNTLIST, min=1)
        validatorNTCHR  = QtPage.IntValidator(self.lineEditNTCHR, min=1)
        validatorNTHIST = QtPage.IntValidator(self.lineEditHisto, min=1)
        validatorFRCHR  = QtPage.DoubleValidator(self.lineEditFRCHR)
        validatorFRHIST = QtPage.DoubleValidator(self.lineEditFRHisto)
        self.lineEditNTLIST.setValidator(validatorNTLIST)
        self.lineEditNTCHR.setValidator(validatorNTCHR)
        self.lineEditHisto.setValidator(validatorNTHIST)
        self.lineEditFRCHR.setValidator(validatorFRCHR)
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

        # Initialisation of the post-pro output frequency

        m = self.mdl.getPostprocessingType()
        if m == 'Frequency_c_x' :
            frchr = self.mdl.getPostprocessingFrequencyTime()
            self.lineEditFRCHR.setText(QString(str(frchr)))
        else :
            ntchr = self.mdl.getPostprocessingFrequency()
            self.lineEditNTCHR.setText(QString(str(ntchr)))
        self.modelPostProcessing.setItem(str_model=m)
        t = self.modelPostProcessing.dicoM2V[m]
        self.slotOutputPostpro(t)

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

        # ICHRVL, ICHRBO, ICHRSY

        if self.mdl.getFluidDomainPostProStatus() == "on":
            self.checkBoxICHRVL.setChecked(True)
        else:
            self.checkBoxICHRVL.setChecked(False)

        if self.mdl.getDomainBoundaryPostProStatus() == "on":
            self.checkBoxICHRBO.setChecked(True)
        else:
            self.checkBoxICHRBO.setChecked(False)

        if ConjugateHeatTransferModel(self.case).getSyrthesCouplingList():
            if self.mdl.getSyrthesBoundaryPostProStatus() == 'on':
                self.checkBoxICHRSY.setChecked(True)
            else:
                self.checkBoxICHRSY.setChecked(False)
        else:
            self.labelICHRSY.hide()
            self.checkBoxICHRSY.hide()

        # values of type of mesh's post processing

        if self.mdl.getTypePostMeshes() == '0':
            if ale == 'on':
                self.modelMeshes.setItem(str_model='10')
            else:
                self.modelMeshes.setItem(str_model='0')
        else:
            self.modelMeshes.setItem(str_model=self.mdl.getTypePostMeshes())


        # values of post processing's format

        fmt = self.mdl.getPostProFormat()
        self.modelFMTCHR.setItem(str_model=fmt)
        line = self.mdl.getPostProOptionsFormat()
        self.__updateOptionsFormat(line)

        # Monitoring points initialisation

        for n in range(self.mdl.getNumberOfMonitoringPoints()):
            name = str(n+1)
            X, Y, Z = self.mdl.getMonitoringPointCoordinates(name)
            self.__insertMonitoringPoint(name, X, Y, Z)


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


    @pyqtSignature("const QString &")
    def slotOutputPostpro(self, text):
        """
        INPUT choice of the output for the Postprocessor (Ensight,...)
        """
        chrono = self.modelPostProcessing.dicoV2M[str(text)]
        log.debug("slotOutputPostpro-> chrono = %s" % chrono)
        self.mdl.setPostprocessingType(chrono)

        if chrono == "At the end":
            ntchr = -1
            self.mdl.setPostprocessingFrequency(ntchr)
            self.lineEditNTCHR.show()
            self.lineEditNTCHR.setText(QString(str(ntchr)))
            self.lineEditNTCHR.setEnabled(False)
            self.lineEditFRCHR.hide()

        elif chrono == "At each step":
            ntchr = 1
            self.mdl.setPostprocessingFrequency(ntchr)
            self.lineEditNTCHR.show()
            self.lineEditNTCHR.setText(QString(str(ntchr)))
            self.lineEditNTCHR.setEnabled(False)
            self.lineEditFRCHR.hide()

        elif chrono == "Frequency_c":
            self.lineEditNTCHR.show()
            self.lineEditNTCHR.setEnabled(True)
            ntchr = self.mdl.getPostprocessingFrequency()
            if ntchr < 1:
                ntchr = 1
                self.mdl.setPostprocessingFrequency(ntchr)
            self.lineEditNTCHR.setText(QString(str(ntchr)))
            self.lineEditFRCHR.hide()

        elif chrono == "Frequency_c_x":
            self.lineEditNTCHR.hide()
            self.lineEditFRCHR.show()
            frchr = self.mdl.getPostprocessingFrequencyTime()
            self.lineEditFRCHR.setText(QString(str(frchr)))


    @pyqtSignature("const QString &")
    def slotPostproFrequency(self, text):
        """
        Input the frequency of the post-processing output
        """
        n, ok = text.toInt()
        if self.sender().validator().state == QValidator.Acceptable:
            log.debug("slotPostproFrequency-> NTCHR = %s" % n)
            self.mdl.setPostprocessingFrequency(n)


    @pyqtSignature("const QString &")
    def slotPostproFrequencyTime(self, text):
        """
        Input the frequency of the post-processing output
        """
        n, ok = text.toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
            log.debug("slotPostproFrequencyTime-> FRCHR = %s" % n)
            self.mdl.setPostprocessingFrequencyTime(n)


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
            self.lineEditHisto.show()
            self.lineEditHisto.setText(QString(str(nthist)))
            self.lineEditHisto.setEnabled(False)
            self.lineEditFRHisto.hide()

        if histo == "At each step":
            nthist = 1
            self.mdl.setMonitoringPointFrequency(nthist)
            self.lineEditHisto.show()
            self.lineEditHisto.setText(QString(str(nthist)))
            self.lineEditHisto.setEnabled(False)
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
    def slotMonitoringPointFrequency(self, text):
        """
        Input the frequency of the monitoring point output
        """
        n, ok = text.toInt()
        if self.sender().validator().state == QValidator.Acceptable:
            log.debug("slotMonitoringPointFrequency-> NTHIST = %s" % n)
            self.mdl.setMonitoringPointFrequency(n)


    @pyqtSignature("const QString &")
    def slotMonitoringPointFrequencyTime(self, text):
        """
        Input the frequency of the monitoring point output
        """
        n, ok = text.toDouble()
        if self.sender().validator().state == QValidator.Acceptable:
            log.debug("slotMonitoringPointFrequencyTime-> FRHIST = %s" % n)
            self.mdl.setMonitoringPointFrequencyTime(n)


    @pyqtSignature("")
    def slotVolumeOuput(self):
        """
        Input value of ICHRVL
        """
        if self.checkBoxICHRVL.isChecked():
            self.mdl.setFluidDomainPostProStatus('on')
        else:
            self.mdl.setFluidDomainPostProStatus('off')


    @pyqtSignature("")
    def slotBoundaryOuput(self):
        """
        Input value of ICHRBO
        """
        if self.checkBoxICHRBO.isChecked():
            self.mdl.setDomainBoundaryPostProStatus('on')
        else:
            self.mdl.setDomainBoundaryPostProStatus('off')

        self.browser.configureTree(self.case)


    @pyqtSignature("")
    def slotSyrthesOutput(self):
        """
        Input value of ICHRSY
        """
        if self.checkBoxICHRSY.isChecked():
            self.mdl.setSyrthesBoundaryPostProStatus('on')
        else:
            self.mdl.setSyrthesBoundaryPostProStatus('off')


    @pyqtSignature("const QString &")
    def slotTypePostMeshes(self, text):
        """
        Input type of post-processing for mesh
        """
        self.mdl.setTypePostMeshes(self.modelMeshes.dicoV2M[str(text)])


    @pyqtSignature("const QString &")
    def slotOutputFormat(self, text):
        """
        Input format of post-processing
        """
        format = self.modelFMTCHR.dicoV2M[str(text)]
        if self.mdl.getPostProFormat() != format:
            self.mdl.setPostProFormat(format)
            line = self.mdl.defaultInitialValues()['postprocessing_options']
            self.mdl.setPostProOptionsFormat(line)
            self.__updateOptionsFormat(line)


    @pyqtSignature("")
    def slotOutputOptions(self):
        """
        Create characters ligne for command of format's options
        """
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
        self.mdl.setPostProOptionsFormat(l)


    def __updateOptionsFormat(self, line):
        """
        Update ligne for command of format's options at each modification of
        post processing format
        """
        list = string.split(line, ',')
        format = self.modelFMTCHR.dicoV2M[str(self.comboBoxFMTCHR.currentText())]
        log.debug("__updateOptionsFormat-> FMTCHR = %s" % format)
        log.debug("__updateOptionsFormat-> OPTCHR = %s" % line)

        # update widgets from the options list

        for opt in list:

            if opt == 'binary' or opt == 'text' :
                self.modelFormat.setItem(str_model=opt)

            if opt == 'discard_polygons' or opt == 'divide_polygons':
                self.modelPolygon.setItem(str_model=opt)

            if opt == 'discard_polyhedra' or opt == 'divide_polyhedra':
                self.modelPolyhedra.setItem(str_model=opt)

            if format == 'EnSight':
                if opt == 'big_endian':
                    self.checkBoxBigEndian.setChecked(True)
##                if opt == 'split_tensors':
##                    self.opt_splittens.set('on')

        if 'discard_polygons' not in list and 'divide_polygons' not in list:
            self.modelPolygon.setItem(str_model="display")
        if 'discard_polyhedra' not in list and 'divide_polyhedra' not in list:
            self.modelPolyhedra.setItem(str_model="display")
        if 'big_endian' not in list:
            self.checkBoxBigEndian.setChecked(False)
##        if 'split_tensors' not in list:
##            self.opt_splittens.set('off')

        # enable and disable options related to the format

        if format != "EnSight":

            if format == "CGNS":
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


    def ale(self):
        """
        """
        ale = MobileMeshModel(self.case).getMethod()
        if ale == 'off':
            self.modelMeshes.disableItem(str_model='10')
            self.modelMeshes.disableItem(str_model='11')
            self.modelMeshes.disableItem(str_model='12')
        return ale


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
