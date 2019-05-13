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
This module defines the main application classes for the Qt GUI.

This module defines the following classes:
- MainView

    @copyright: 1998-2017 EDF S.A., France
    @author: U{EDF<mailto:saturne-support@edf.fr>}
    @license: GNU GPL v2, see COPYING for details.
"""

#-------------------------------------------------------------------------------
# Standard modules
#-------------------------------------------------------------------------------

import os, sys, shutil, signal, logging
import subprocess, platform

try:
    import ConfigParser  # Python2
    configparser = ConfigParser
except Exception:
    import configparser  # Python3

#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------

from code_saturne.Base.QtCore    import *
from code_saturne.Base.QtGui     import *
from code_saturne.Base.QtWidgets import *

#-------------------------------------------------------------------------------
# Application modules
#-------------------------------------------------------------------------------

import cs_info
from cs_exec_environment import \
    separate_args, update_command_single_value, assemble_args, enquote_arg
import cs_runcase

try:
    from code_saturne.trackcvg.MainForm import Ui_MainForm
except:
    sys.path.insert(1, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "Base"))
    from code_saturne.trackcvg.MainForm import Ui_MainForm

try:
    import code_saturne.trackcvg
except:
    sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from code_saturne.Base.QtPage import getexistingdirectory
from code_saturne.Base.QtPage import DoubleValidator, from_qvariant

import numpy
import matplotlib
import matplotlib.pyplot
if QT_API == 'PYQT4':
    import matplotlib.backends.backend_qt4agg
    from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
    from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
else:
    import matplotlib.backends.backend_qt5agg
    from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
    from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure

import xml
from xml.dom.minidom import parse, Document

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("MainView")

#-------------------------------------------------------------------------------
# item class
#-------------------------------------------------------------------------------

class item_class(object):
    '''
    custom data object
    '''
    def __init__(self, idx, name, status, subplot_id):
        nm, ext = os.path.splitext(name)
        self.index      = idx
        self.name       = nm
        self.status     = status
        self.subplot_id = subplot_id

    def __repr__(self):
        return "case : %s // status %s"\
               % (self.name, self.status)


#-------------------------------------------------------------------------------
# Treeitem class
#-------------------------------------------------------------------------------

class TreeItem(object):
    '''
    a python object used to return row/column data, and keep note of
    it's parents and/or children
    '''
    def __init__(self, item, header, parentItem):
        self.item = item
        self.parentItem = parentItem
        self.header = header
        self.childItems = []


    def appendChild(self, item):
        self.childItems.append(item)


    def child(self, row):
        return self.childItems[row]


    def childCount(self):
        return len(self.childItems)


    def columnCount(self):
        return 3


    def data(self, column, role):
        if self.item == None:
            if column == 0:
                self.header
            else:
                return None
        else:
            if column == 0 and role == Qt.DisplayRole:
                return self.item.name
            elif column == 1 and role == Qt.CheckStateRole:
                value = self.item.status
                if value == 'on':
                    return Qt.Checked
                elif value == 'onoff':
                    return Qt.PartiallyChecked
                else:
                    return Qt.Unchecked
            elif column == 2 and role == Qt.DisplayRole:
                return self.item.subplot_id
        return None


    def parent(self):
        return self.parentItem


    def row(self):
        if self.parentItem:
            return self.parentItem.childItems.index(self)
        return 0


#-------------------------------------------------------------------------------
# StandarItemModelOutput class
#-------------------------------------------------------------------------------

class CaseStandardItemModel(QAbstractItemModel):

    def __init__(self, parent, lst, lstFileProbes):
        """
        """
        QAbstractItemModel.__init__(self)

        self.lst           = lst
        self.lstFileProbes = lstFileProbes

        self.noderoot = {}

        self.rootItem = TreeItem(None, "ALL", None)
        self.parents = {0 : self.rootItem}

        self.populateModel()


    def columnCount(self, parent = None):
        if parent and parent.isValid():
            return parent.internalPointer().columnCount()
        else:
            return 3


    def data(self, index, role):
        if not index.isValid():
            return None

        item = index.internalPointer()

        # ToolTips
        if role == Qt.ToolTipRole:
            return None

        # StatusTips
        if role == Qt.StatusTipRole:
            if index.column() == 0:
                return self.tr("File name")
            elif index.column() == 1:
                return self.tr("Status")
            elif index.column() == 2:
                return self.tr("Subplot")

        # Display
        if role == Qt.DisplayRole:
            return item.data(index.column(), role)
        elif role == Qt.CheckStateRole:
            return item.data(index.column(), role)

        return None


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled

        itm = index.internalPointer()

        if itm in self.noderoot.values():
            # traitement des categories
            if index.column() == 0:
                return Qt.ItemIsEnabled | Qt.ItemIsSelectable
            elif index.column() == 1:
                return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsUserCheckable | Qt.ItemIsTristate
            elif index.column() == 2:
                return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable
        else:
            if index.column() == 0:
                return Qt.ItemIsEnabled | Qt.ItemIsSelectable
            elif index.column() == 1:
                return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsUserCheckable
            elif index.column() == 2:
                return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable


    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            if section == 0:
                return self.tr("File name")
            elif section == 1:
                return self.tr("Status")
            elif section == 2:
                return self.tr("Subplot")
        return None


    def index(self, row, column, parent = QModelIndex()):
        if not self.hasIndex(row, column, parent):
            return QModelIndex()

        if not parent.isValid():
            parentItem = self.rootItem
        else:
            parentItem = parent.internalPointer()

        try:
            childItem = parentItem.child(row)
        except:
            childItem = None

        if childItem:
            return self.createIndex(row, column, childItem)
        else:
            return QModelIndex()


    def parent(self, index):
        if not index.isValid():
            return QModelIndex()

        childItem = index.internalPointer()
        if not childItem:
            return QModelIndex()

        parentItem = childItem.parent()

        if parentItem == self.rootItem:
            return QModelIndex()

        return self.createIndex(parentItem.row(), 0, parentItem)


    def rowCount(self, parent=QModelIndex()):
        if not parent.isValid():
            p_Item = self.rootItem
        else:
            p_Item = parent.internalPointer()
        return p_Item.childCount()


    def populateModel(self):
        idx = 0
        for (name, name_long, status, subplot_id, probes_number) in self.lst:
            item = item_class(idx, name, status, subplot_id)
            newparent = TreeItem(item, name, self.rootItem)
            self.rootItem.appendChild(newparent)
            self.noderoot[name] = newparent
            idx = idx + 1

        # create subnodes
        for (name, name_long, status, subplot_id, probes_number) in self.lst:
            ll = []
            for itm in self.lstFileProbes[name]:  # first column is time or iteration
                parentItem = self.noderoot[name]
                nameItem = itm.name
                newItem = TreeItem(itm, nameItem, parentItem)
                parentItem.appendChild(newItem)


    def setData(self, index, value, role=None):
        item = index.internalPointer()

        if index.column() == 1:
            v = from_qvariant(value, int)
            if v == Qt.Checked:
                item.item.status = "on"
            else:
                item.item.status = "off"
            if item in self.noderoot.values():
                (name, name_long, status,
                 subplot_id, probes_number) = self.lst[index.row()]
                self.lst[index.row()] = (name, name_long, item.item.status,
                                         subplot_id, probes_number)
                for it in self.lstFileProbes[name]:
                    it.status = item.item.status
            else:
                size = len(item.parentItem.childItems)
                active = 0
                for itm in item.parentItem.childItems:
                    if itm.item.status == "on":
                        active = active + 1
                if active == 0:
                    item.parentItem.item.status = "off"
                elif active == size:
                    item.parentItem.item.status = "on"
                else:
                    item.parentItem.item.status = "onoff"
                # update self.lst
                idx = 0
                for (name, name_long, status, subplot_id, probes_number) in self.lst:
                    nm, ext = os.path.splitext(name)
                    if nm == item.parentItem.item.name:
                        self.lst[idx] = (name, name_long,
                                         item.parentItem.item.status,
                                         subplot_id, probes_number)
                    idx = idx + 1

        elif index.column() == 2:
            v = from_qvariant(value, int)
            item.item.subplot_id = v
            if item in self.noderoot.values():
                (name, name_long, status,
                 subplot_id, probes_number) = self.lst[index.row()]
                self.lst[index.row()] = (name, name_long, status,
                                         item.item.subplot_id, probes_number)
                for it in self.lstFileProbes[name]:
                    it.subplot_id = item.item.subplot_id
            else:
                size = len(item.parentItem.childItems)
                subplotid = item.parentItem.childItems[0]
                for itm in item.parentItem.childItems:
                    if itm.item.subplot_id != subplotid:
                        subplotid = -1
                # update self.lst
                item.parentItem.item.subplot_id = subplotid
                idx = 0
                for (name, name_long, status, subplot_id, probes_number) in self.lst:
                    nm, ext = os.path.splitext(name)
                    if nm == item.parentItem.item.name:
                        self.lst[idx] = (name, name_long, status,
                                         subplotid, probes_number)
                    idx = idx + 1

        self.dataChanged.emit(QModelIndex(), QModelIndex())

        return True


#-------------------------------------------------------------------------------
# manage figures
#-------------------------------------------------------------------------------

class MyMplCanvas(FigureCanvas):
    """Ultimately, this is a QWidget (as well as a FigureCanvasAgg, etc.)."""

    def __init__(self, parent=None, subplotNb=1, width=5, height=4, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        super(MyMplCanvas,self).__init__(self.fig)
        self.yAxe = numpy.array([0])
        self.xAxe = numpy.array([0])
        self.axes = []

        if subplotNb == 1:
            self.axes.append(self.fig.add_subplot(111))
        elif subplotNb == 2:
            self.axes.append(self.fig.add_subplot(211))
            self.axes.append(self.fig.add_subplot(212))
        elif subplotNb == 3:
            self.axes.append(self.fig.add_subplot(221))
            self.axes.append(self.fig.add_subplot(222))
            self.axes.append(self.fig.add_subplot(223))
        elif subplotNb == 4:
            self.axes.append(self.fig.add_subplot(221))
            self.axes.append(self.fig.add_subplot(222))
            self.axes.append(self.fig.add_subplot(223))
            self.axes.append(self.fig.add_subplot(224))

        #self.axes.autoscale(False)
        # We want the axes cleared every time plot() is called
        #self.axes.hold(False)

        self.compute_initial_figure()

        #
        FigureCanvas.__init__(self, self.fig)
        self.setParent(parent)

        FigureCanvas.setSizePolicy(self,
                                   QSizePolicy.Expanding,
                                   QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)


    def compute_initial_figure(self):
        pass


    def update_figure(self, name, data, nb_probes, lstProbes):
        self.xAxe = data[0].tolist()
        for j in range(nb_probes - 1):
            if (lstProbes[j].status == "on"):
                self.yAxe = data[j + 1].tolist()

                lbl = name + "_s" + str(j)

                self.axes[lstProbes[j].subplot_id - 1].plot(self.xAxe, self.yAxe,
                                                            label = lbl)
                self.axes[lstProbes[j].subplot_id - 1].legend(loc="upper left",
                                                              bbox_to_anchor=(1.02,1.0),
                                                              borderaxespad=0.0,
                                                              ncol=1,
                                                              fancybox=True,
                                                              shadow=True,
                                                              prop={'size':'medium',
                                                                    'style': 'italic'})


    def update_figure_listing(self, name, data, nb_probes, lstProbes):
        self.xAxe = data[0].tolist()
        for j in range(nb_probes - 1):
            if (lstProbes[j].status == "on"):
                self.yAxe = data[j + 1].tolist()

                lbl = "t res. " + name[j]

                self.axes[lstProbes[j].subplot_id - 1].plot(self.xAxe, self.yAxe,
                                                            label = lbl)
                self.axes[lstProbes[j].subplot_id - 1].legend(loc="upper left",
                                                              bbox_to_anchor=(1.02,1.0),
                                                              borderaxespad=0.0,
                                                              ncol=1,
                                                              fancybox=True,
                                                              shadow=True,
                                                              prop={'size':'medium',
                                                                    'style': 'italic'})



    def drawFigure(self):
        for it in range(len(self.axes)):
            self.axes[it].grid(True)
            self.axes[it].set_xlabel("time (s)")
        self.axes[0].set_yscale('log')

        self.fig.canvas.draw()


    def clear(self):
        for plt in range(len(self.axes)):
            self.axes[plt].clear()


    def setSubplotNumber(self, subplotNb):
        self.fig.clear()
        for it in range(len(self.axes)):
            self.axes.remove(self.axes[0])
        if subplotNb == 1:
            self.axes.append(self.fig.add_subplot(111))
        elif subplotNb == 2:
            self.axes.append(self.fig.add_subplot(211))
            self.axes.append(self.fig.add_subplot(212))
        elif subplotNb == 3:
            self.axes.append(self.fig.add_subplot(221))
            self.axes.append(self.fig.add_subplot(222))
            self.axes.append(self.fig.add_subplot(223))
        elif subplotNb == 4:
            self.axes.append(self.fig.add_subplot(221))
            self.axes.append(self.fig.add_subplot(222))
            self.axes.append(self.fig.add_subplot(223))
            self.axes.append(self.fig.add_subplot(224))


#-------------------------------------------------------------------------------
# Base Main Window
#-------------------------------------------------------------------------------

class MainView(object):
    """
    Abstract class
    """
    NextId = 1
    Instances = set()

    def __new__(cls, cmd_package = None, cmd_case = ""):
        """
        Factory
        """
        return MainViewSaturne.__new__(MainViewSaturne, cmd_package, cmd_case)


    @staticmethod
    def updateInstances(qobj):
        """
        Overwrites the Instances set with a set that contains only those
        window instances that are still alive.
        """
        MainView.Instances = set([window for window \
                in MainView.Instances if isAlive(window)])


    def ui_initialize(self):
        self.setAttribute(Qt.WA_DeleteOnClose)
        MainView.Instances.add(self)

        iconpath = os.path.split(os.path.dirname(os.path.abspath(__file__)))[0]
        iconpath = os.path.join(iconpath, "Base", "MONO-bulle-HD.png")
        icon = QIcon(QPixmap(iconpath))
        self.setWindowIcon(icon)

        self.setWindowTitle(self.package.code_name + " TRACKING CONVERGENCE" + " - " + self.package.version)

        # Validator
        validator = DoubleValidator(self.lineEditTime, min=0.0)
        self.lineEditTime.setValidator(validator)

        # connections
        self.fileCloseAction.triggered.connect(self.caseClose)
        self.fileQuitAction.triggered.connect(self.fileQuit)
        self.actionSave_state.triggered.connect(self.SaveState)
        self.actionLoad_state.triggered.connect(self.LoadState)
        self.actionSave_state.setEnabled(False)
        self.actionLoad_state.setEnabled(False)

        self.displayAboutAction.triggered.connect(self.displayAbout)
        self.backgroundColorAction.triggered.connect(self.setColor)
        self.actionFont.triggered.connect(self.setFontSize)
        self.RestoreStyleDefaults.triggered.connect(self.restoreStyleDefaults)

        self.displayLicenceAction.triggered.connect(self.displayLicence)

        self.toolButtonDir.clicked.connect(self.slotOpenCase)
        self.lineEditTime.textChanged[str].connect(self.slotRefreshTime)
        self.pushButtonRefresh.clicked.connect(self.slotRefresh)

        # connection for page layout
        self.destroyed.connect(MainView.updateInstances)

        # Ctrl+C signal handler (allow to shutdown the GUI with Ctrl+C)

        signal.signal(signal.SIGINT, signal.SIG_DFL)

        self.resize(800, 700)

        # restore system settings

        settings = QSettings()

        try:
            self.restoreGeometry(settings.value("MainWindow/Geometry", QByteArray()))
            self.restoreState(settings.value("MainWindow/State", QByteArray()))
        except:
            self.restoreGeometry(settings.value("MainWindow/Geometry").toByteArray())
            self.restoreState(settings.value("MainWindow/State").toByteArray())

        app = QCoreApplication.instance()

        self.palette_default = None
        self.font_default = None

        if settings.contains("MainWindow/Color"):
            color = settings.value("MainWindow/Color",
                                   self.palette().color(QPalette.Window).name())
            color = QColor(color)
            if color.isValid():
                if not self.palette_default:
                    self.palette_default = QPalette().resolve(self.palette())
                self.setPalette(QPalette(color))
                app.setPalette(QPalette(color))

        if settings.contains("MainWindow/Font"):
            f = settings.value("MainWindow/Font", str(self.font()))
            if f:
                if not self.font_default:
                    self.font_default = self.font()
                font = QFont()
                if (font.fromString(from_qvariant(f, to_text_string))):
                    self.setFont(font)
                    app.setFont(font)

        # Init
        if self.caseName == None:
            self.toolButtonDir.setStyleSheet("background-color: red")
        else:
            self.openCase(self.caseName)
        self.lineEditCase.setEnabled(False)
        self.lineEditTime.setText(str(self.timeRefresh))

        # treeViewDirectory
        self.modelCases = CaseStandardItemModel(self.parent,
                                                self.fileList,
                                                self.listFileProbes)
        self.treeViewDirectory.setModel(self.modelCases)
        self.treeViewDirectory.setAlternatingRowColors(True)
        self.treeViewDirectory.setSelectionBehavior(QAbstractItemView.SelectItems)
        self.treeViewDirectory.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.treeViewDirectory.setEditTriggers(QAbstractItemView.DoubleClicked)
        self.treeViewDirectory.expandAll()
        self.treeViewDirectory.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.treeViewDirectory.setDragEnabled(False)
        self.treeViewDirectory.resizeColumnToContents(0)
        self.treeViewDirectory.resizeColumnToContents(1)
        self.treeViewDirectory.resizeColumnToContents(2)
        self.treeViewDirectory.resizeColumnToContents(3)
        self.modelCases.dataChanged.connect(self.treeViewChanged)

        self.timer.timeout.connect(self.updateView)

        self.spinBox.valueChanged[int].connect(self.updateSubplotNumber)

        # gestion des figures
        l = QVBoxLayout(self.widget)
        self.dc = MyMplCanvas(self.widget, subplotNb=self.subplotNumber,
                              width=5, height=4, dpi=50)
        l.addWidget(self.dc)

        # this is the Navigation widget
        # it takes the Canvas widget and a parent
        self.toolbar = NavigationToolbar(self.dc, self.widget)
        l.addWidget(self.toolbar)

        self.statusbar.setSizeGripEnabled(False)
        self.statusbar.showMessage(self.tr("Ready"), 5000)
        self.updateView()


    def caseClose(self):
        """
        public slot

        try to quit all the current MainWindow
        """
        if self.caseName != None:
            title = self.tr("Close Case")
            msg   = self.tr("Save current state?")
            reply = QMessageBox.question(self, title, msg,
                                         QMessageBox.Yes | QMessageBox.No)
            if reply == QMessageBox.Yes:
                self.SaveState()

        self.toolButtonDir.setStyleSheet("background-color: red")
        self.lineEditCase.setText("")
        self.caseName = None
        self.fileList = []
        self.listingVariable = []
        self.listFileProbes = {}
        self.modelCases = CaseStandardItemModel(self.parent, [], [])
        self.treeViewDirectory.setModel(self.modelCases)
        self.modelCases.dataChanged.connect(self.treeViewChanged)
        self.updateView()


    def closeEvent(self, event):
        """
        public slot

        try to quit all the current MainWindow
        """
        if self.caseName != None:
            self.caseClose()

        settings = QSettings()

        settings.setValue("MainWindow/Geometry",
                          self.saveGeometry())
        settings.setValue("MainWindow/State",
                          self.saveState())

        event.accept()


    def fileQuit(self):
        """
        Public slot.

        try to quit all window
        """
        QApplication.closeAllWindows()


    def displayManual(self, pkg, manual, reader = None):
        """
        private method

        open a manual
        """
        argv_info = ['--guide']
        argv_info.append(manual)
        cs_info.main(argv_info, pkg)


    def displayAbout(self):
        """
        public slot

        the About dialog window shows:
         - title
         - version
         - contact
        """
        msg = self.package.code_name + "\n"                 +\
              "version " + self.package.version + "\n\n"    +\
              "For information about this application "     +\
              "please contact:\n\n"                         +\
              self.package.bugreport + "\n\n"               +\
              "Please visit our site:\n"                    +\
              self.package.url
        QMessageBox.about(self, self.package.name + ' convergence plot', msg)


    def displayLicence(self):
        """
        public slot

        GNU GPL license dialog window
        """
        QMessageBox.about(self, self.package.code_name + ' convergence plot',
                          cs_info.licence_text)


    def displayConfig(self):
        """
        public slot

        configuration information window
        """
        QMessageBox.about(self, self.package.code_name + ' convergence plot',
                          "see config.py")


    def setColor(self):
        """
        public slot

        choose GUI color
        """
        c = self.palette().color(QPalette.Window)
        color = QColorDialog.getColor(c, self)
        if color.isValid():
            app = QCoreApplication.instance()
            if not self.palette_default:
                self.palette_default = QPalette().resolve(self.palette())
            app.setPalette(QPalette(color))
            settings = QSettings()
            settings.setValue("MainWindow/Color",
                              self.palette().color(QPalette.Window).name())


    def setFontSize(self):
        """
        public slot

        choose GUI font
        """
        font, ok = QFontDialog.getFont(self)
        log.debug("setFont -> %s" % ok)
        if ok:
            if not self.font_default:
                self.font_default = self.font()
            self.setFont(font)
            app = QCoreApplication.instance()
            app.setFont(font)
            settings = QSettings()
            settings.setValue("MainWindow/Font",
                              self.font().toString())


    def restoreStyleDefaults(self):
        """
        public slot

        Restore default style.
        """

        reply = QMessageBox.question(self, "Restore defaults",
                                     "Restore default color and font ?",
                                     QMessageBox.Yes | QMessageBox.No)
        if reply == QMessageBox.Yes:
            app = QCoreApplication.instance()
            if self.palette_default:
                app.setPalette(self.palette_default)
            if self.font_default:
                print(self.font_default)
                print(self.font())
                self.setFont(self.font_default)
                app.setFont(self.font_default)
            settings = QSettings()
            settings.remove("MainWindow/Color")
            settings.remove("MainWindow/Font")


    def loadDirectoryContent(self):
        """
        Load directory content in treeView
        """
        # self.fileList : nom_court, nom_long, status, plotId, probes_number
        for fl in os.listdir(self.caseName):
            rep = os.path.abspath(os.path.join(self.caseName, fl))
            if os.path.isdir(rep):
                for ffl in os.listdir(rep):
                    base, ext = os.path.splitext(ffl)
                    if ext in ['.dat', '.csv']:
                        # read number of probes
                        if ext == ".csv" and (base.find("_coords") == -1):
                            size = self.ReadCsvFileHeader(os.path.abspath(os.path.join(rep, ffl)))
                            self.fileList.append([ffl, os.path.abspath(os.path.join(rep, ffl)), "off", 2, size])
                            ll = []
                            for idx in range(size - 1):
                                nameItem = "probe_" + str(idx)
                                item = item_class(idx, nameItem, "off", 2)
                                ll.append(item)
                            self.listFileProbes[ffl] = ll
                        elif ext == ".dat" and (base.find("_coords") == -1):
                            size = self.ReadDatFileHeader(os.path.abspath(os.path.join(rep, ffl)))
                            self.fileList.append([ffl, os.path.abspath(os.path.join(rep, ffl)), "off", 2, size])
                            ll = []
                            for idx in range(size - 1):
                                nameItem = "probe_" + str(idx)
                                item = item_class(idx, nameItem, "off", 2)
                                ll.append(item)
                            self.listFileProbes[ffl] = ll
            else:
                if fl == 'residuals.csv':
                    self.listingVariable = self.readResidualsVariableListCSV(rep)
                    self.fileList.append([fl, rep, "on", 1, len(self.listingVariable) + 1])
                    # read variable list for Time residual
                    idx = 0
                    ll = []
                    for var in self.listingVariable:
                        item = item_class(idx, var, "on", 1)
                        idx = idx + 1
                        ll.append(item)
                    self.listFileProbes[fl] = ll
                elif fl == 'residuals.dat':
                    self.listingVariable = self.readResidualsVariableListDAT(rep)
                    self.fileList.append([fl, rep, "on", 1, len(self.listingVariable) + 1])
                    # read variable list for Time residual
                    idx = 0
                    ll = []
                    for var in self.listingVariable:
                        item = item_class(idx, var, "on", 1)
                        idx = idx + 1
                        ll.append(item)
                    self.listFileProbes[fl] = ll
        self.LoadState()

        self.modelCases = CaseStandardItemModel(self.parent, self.fileList, self.listFileProbes)
        self.treeViewDirectory.setModel(self.modelCases)
        self.treeViewDirectory.expandAll()
        self.treeViewDirectory.resizeColumnToContents(0)
        self.treeViewDirectory.resizeColumnToContents(1)
        self.treeViewDirectory.resizeColumnToContents(2)
        self.treeViewDirectory.resizeColumnToContents(3)
        self.modelCases.dataChanged.connect(self.treeViewChanged)


    def treeViewChanged(self, topLeft, bottomRight):
        """
        """
        self.updateView()


    def openCase(self, dirName):

        self.caseName = str(dirName)
        self.toolButtonDir.setStyleSheet("background-color: green")
        self.lineEditCase.setText(self.caseName)
        self.loadDirectoryContent()

        if len(self.fileList) > 0:
            self.actionSave_state.setEnabled(True)
        else:
            self.actionSave_state.setEnabled(False)

        if os.path.exists(os.path.join(self.caseName, '.trackcvg.state')):
            self.actionLoad_state.setEnabled(True)
        else:
            self.actionLoad_state.setEnabled(False)


    def slotOpenCase(self):

        # Save previous case before switching

        if self.caseName != None:
            self.caseClose()

        # Load new case

        title = self.tr("Choose a result directory")

        path = os.getcwd()
        dataPath = os.path.join(path, "..", "RESU")
        if os.path.isdir(dataPath): path = dataPath

        dirName = getexistingdirectory(self, title,
                                       path,
                                       QFileDialog.ShowDirsOnly |
                                       QFileDialog.DontResolveSymlinks)
        if not dirName:
            self.caseName = None
            return

        self.openCase(dirName)


    def slotRefreshTime(self, text):
        """
        Private slot.
        """
        if self.lineEditTime.validator().state == QValidator.Acceptable:
            time = from_qvariant(text, float)
            self.timeRefresh = time
            if self.timeRefresh != 0.:
                self.timer.start(self.timeRefresh * 1000)
            else:
                self.timer.start(1.e8 * 1000)


    def updateSubplotNumber(self, v):
        """
        Increment, decrement number of subplot
        """
        old_num = self.subplotNumber
        self.subplotNumber = int(self.spinBox.text())
        self.dc.setSubplotNumber(self.subplotNumber)
        if self.subplotNumber < old_num:
            for it in self.listFileProbes:
                for itt in self.listFileProbes[it]:
                    if itt.subplot_id > self.subplotNumber:
                        itt.subplot_id = 1
                        itt.status = "off"
        self.updateView()
        # TODO manque update au niveau fichier!!!!!! et le traitement de NEPTUNE


    def updateView(self):
        """
        """
        self.dc.clear()
        for (name, fle, status, subplot_id, probes_number) in self.fileList:
            if status == "on" or status == "onoff":
                base, ext = os.path.splitext(fle)
                if name == 'residuals.csv':
                    data = self.ReadCsvFile(fle, probes_number)
                    if status == "on" or status == "onoff":
                        self.dc.update_figure_listing(self.listingVariable,
                                                      data,
                                                      probes_number,
                                                      self.listFileProbes[name])
                elif name == 'residuals.dat':
                    data = self.ReadDatFile(fle, probes_number)
                    if status == "on" or status == "onoff":
                        self.dc.update_figure_listing(self.listingVariable,
                                                      data, probes_number,
                                                      self.listFileProbes[name])
                elif ext == ".csv":
                    data = self.ReadCsvFile(fle, probes_number)
                    nm, ext = os.path.splitext(name)
                    nm = nm[7:]
                    self.dc.update_figure(nm, data, probes_number,
                                          self.listFileProbes[name])
                elif ext == ".dat":
                    data = self.ReadDatFile(fle, probes_number)
                    nm, ext = os.path.splitext(name)
                    nm = nm[7:]
                    self.dc.update_figure(nm, data, probes_number,
                                          self.listFileProbes[name])
        self.dc.drawFigure()


    def ReadCsvFile(self, name, probes_number):
        """
        """
        data = []
        comp = 0
        ficIn= open(name, 'r')
        typ = 0 # O time / 1 iteration
        for line in ficIn.readlines():
            e_id = line.find('\n')
            if e_id < 0:
                continue
            if comp > 0:
                content = line[:e_id].split(',')

                if comp == 1:
                    if type(content[0]) == int:
                        typ = 1

                if typ == 0:
                    for el in content:
                        el = el.lstrip()
                        el = el.rstrip()
                        data.append(el)
                else:
                    data.append(float(content[0]))
                    for el in content[1:]:
                        el = el.lstrip()
                        el = el.rstrip()
                        data.append(el)
            comp = comp + 1
        ficIn.close()
        A = numpy.array(data)
        n = A.shape[0]
        r = n % probes_number
        if r > 0:
            n -= r
            A = numpy.resize(A, n)
        B = A.reshape(-1, probes_number)
        C = B.transpose()

        return C


    def ReadCsvFileHeader(self, name):
        """
        """
        size = 0
        ficIn= open(name, 'r')
        # first line x, y, z
        line = ficIn.readline()
        line = line.strip('\n')
        content = line.split(',')
        ficIn.close()

        return len(content)


    def ReadDatFile(self, name, probes_number):
        """
        """
        data = []
        ficIn= open(name, 'r')
        comp = 0
        typ = 0 # O time / 1 iteration
        for line in ficIn.readlines():
            if not line.startswith("#"):
                comp = comp + 1
                e_id = line.find('\n')
                if e_id < 0:
                    continue
                line = line[:e_id].lstrip()
                content = line.split()

                if comp == 1:
                    if type(content[0]) == int:
                        typ = 1

                if typ == 0:
                    for el in content:
                        data.append(el)
                else:
                    data.append(float(content[0]))
                    for el in content[1:]:
                        data.append(el)
        ficIn.close()
        A = numpy.array(data)
        n = A.shape[0]
        r = n % probes_number
        if r > 0:
            n -= r
            A = numpy.resize(A, n)
        B = A.reshape(-1, probes_number)
        C = B.transpose()

        return C


    def ReadDatFileHeader(self, name):
        """
        """
        size = 0
        ficIn= open(name, 'r')
        for line in ficIn.readlines():
            if not line.startswith("#"):
                line = line.lstrip()
                content = line.split()
                lst = []
                size = len(content)
        ficIn.close()

        return size


    def readResidualsVariableListCSV(self, name):
        """
        """
        ficIn= open(name, 'r')
        lst = []
        line  = ficIn.readline()
        line = line.strip('\n')
        line = line.lstrip()
        content = line.split(',')
        for el in content[1:]:
            lst.append(el)
        ficIn.close()
        return lst


    def readResidualsVariableListDAT(self, name):
        """
        """
        ficIn= open(name, 'r')
        lst = []
        line  = ficIn.readline()
        line = line.lstrip()
        content = line.split()
        for el in content[1:]:
            lst.append(el)
        ficIn.close()
        return lst


    def SaveState(self):
        """
        """
        newdoc = Document()

        root = newdoc.createElement('root')
        newdoc.appendChild(root)

        newnode = newdoc.createElement('timeRefresh')
        newnode.setAttribute("value", str(self.timeRefresh))
        root.appendChild(newnode)
        newnode = newdoc.createElement('subplotNumber')
        newnode.setAttribute("value", str(self.subplotNumber))
        root.appendChild(newnode)

        for (name, fle, status, subplot_id, probes_number) in self.fileList:
            newnode = newdoc.createElement('file')
            newnode.setAttribute("name", name)
            root.appendChild(newnode)

            for itt in self.listFileProbes[name]:
                node = newdoc.createElement('probe')
                node.setAttribute("id", str(itt.index))
                node.setAttribute("name", itt.name)
                node.setAttribute("status", itt.status)
                node.setAttribute("subplot_id", str(itt.subplot_id))
                newnode.appendChild(node)

        name = os.path.join(self.caseName, '.trackcvg.state')
        ficIn= open(name, 'w')

        newdoc.writexml(ficIn,
                        indent="  ",
                        addindent="  ",
                        newl='\n')

        newdoc.unlink()
        ficIn.close()


    def LoadState(self):
        """
        """
        name = os.path.join(self.caseName, '.trackcvg.state')

        if os.path.exists(name):
            dom = parse(name)
            self.timeRefresh = float(dom.getElementsByTagName('timeRefresh')[0].getAttribute('value'))
            self.subplotNumber = int(dom.getElementsByTagName('subplotNumber')[0].getAttribute('value'))
            for (name, fle, status, subplot_id, probes_number) in self.fileList:
                lst = dom.getElementsByTagName('file')
                fl = None
                for node in lst:
                    if node.getAttribute('name') == name:
                        fl = node
                if fl != None:
                    for itt in self.listFileProbes[name]:
                        nn = fl.getElementsByTagName('probe')[itt.index]
                        itt.name       = nn.getAttribute('name')
                        itt.status     = nn.getAttribute('status')
                        itt.subplot_id = int(nn.getAttribute('subplot_id'))
            self.lineEditTime.setText(str(self.timeRefresh))
            self.timer.start(self.timeRefresh * 1000)
            self.spinBox.setValue(int(self.subplotNumber))


    def slotRefresh(self):
        """
        """
        if not self.caseName:
            return

        name = os.path.join(self.caseName, 'control_file')
        ficIn= open(name, 'w')
        ficIn.write('flush\n')
        ficIn.close()


    def tr(self, text):
        """
        private method

        translation

        @param text: text to translate
        @return: translated text
        """
        return text

#-------------------------------------------------------------------------------
# Main Window for Code_Saturne
#-------------------------------------------------------------------------------

class MainViewSaturne(QMainWindow, Ui_MainForm, MainView):

    def __new__(cls, cmd_package = None, cmd_case = ""):
        return super(MainViewSaturne, cls). __new__(cls, cmd_package, cmd_case)


    def __init__(self,
                 cmd_package      = None,
                 cmd_case         = ""):
        """
        Initializes a Main Window for a new document:
          1. finish the Main Window layout
          2. connection betwenn signal and slot
          3. Ctrl+C signal handler
          4. create some instance variables
          5. restore system settings

        @type cmd_case:
        @param cmd_case:
        """
        QMainWindow.__init__(self)
        Ui_MainForm.__init__(self)

        self.setupUi(self)

        # create some instance variables

        self.package     = cmd_package
        if cmd_case != None:
            self.caseName = os.path.abspath(cmd_case)
        else:
            self.caseName = cmd_case
        self.timeRefresh = 10.
        self.subplotNumber = 2
        self.fileList = []
        self.listingVariable = []
        self.listFileProbes = {}
        self.timer = QTimer()
        self.timer.start(self.timeRefresh * 1000)

        self.ui_initialize()

        self.displayCSManualAction.triggered.connect(self.displayCSManual)
        self.displayCSTutorialAction.triggered.connect(self.displayCSTutorial)
        self.displayCSTheoryAction.triggered.connect(self.displayCSTheory)
        self.displayCSSmgrAction.triggered.connect(self.displayCSSmgr)
        self.displayCSRefcardAction.triggered.connect(self.displayCSRefcard)
        self.displayCSDoxygenAction.triggered.connect(self.displayCSDoxygen)

        docdir = self.package.get_dir('docdir')
        if os.path.isdir(docdir):
            liste = os.listdir(docdir)
        else:
            liste = []

        if 'user.pdf' not in liste:
            self.displayCSManualAction.setEnabled(False)
        if 'theory.pdf' not in liste:
            self.displayCSTheoryAction.setEnabled(False)
        if 'studymanager.pdf' not in liste:
            self.displayCSSmgrAction.setEnabled(False)
        if 'refcard.pdf' not in liste:
            self.displayCSRefcardAction.setEnabled(False)
        if 'doxygen' not in liste:
            self.displayCSDoxygenAction.setEnabled(False)
        self.displayNCManualAction.setVisible(False)


    def displayCSManual(self):
        """
        public slot

        open the user manual
        """
        self.displayManual(self.package, 'user')


    def displayCSTutorial(self):
        """
        public slot

        open the tutorial for Code_Saturne
        """
        msg = "See " + self.package.url + " web site for tutorials."
        QMessageBox.about(self, self.package.name + ' convergence plot', msg)


    def displayCSTheory(self):
        """
        public slot

        open the theory and programmer's guide
        """
        self.displayManual(self.package, 'theory')

    def displayCSSmgr(self):
        """
        public slot

        open the studymanager guide
        """
        self.displayManual(self.package, 'studymanager')

    def displayCSRefcard(self):
        """
        public slot

        open the quick reference card for Code_Saturne
        """
        self.displayManual(self.package, 'refcard')


    def displayCSDoxygen(self):
        """
        public slot

        open the quick doxygen for Code_Saturne
        """
        self.displayManual(self.package, 'Doxygen')


#-------------------------------------------------------------------------------

def isAlive(qobj):
    """
    return True if the object qobj exist

    @param qobj: the name of the attribute
    @return: C{True} or C{False}
    """
    import sip
    try:
        sip.unwrapinstance(qobj)
    except RuntimeError:
        return False
    return True

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
