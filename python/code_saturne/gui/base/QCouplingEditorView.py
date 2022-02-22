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
This module defines the view for coupling parameters handling.
"""
#-------------------------------------------------------------------------------
# Standard modules
#-------------------------------------------------------------------------------

import sys, os
import logging
from code_saturne.gui.base.QtCore    import *
from code_saturne.gui.base.QtGui     import *
from code_saturne.gui.base.QtWidgets import *

# Check if QString exists
has_qstring = True
try:
    from code_saturne.gui.base.QtCore import QString
    _fromUtf8 = QString.fromUtf8
except ImportError:
    has_qstring = False
    def _fromUtf8(s):
        return s

    def QString(s):
        return s

try:
    # PyQt5
    from code_saturne.gui.base.QtWidgets import QMainWindow, QMessageBox, \
        QAction, QFileDialog, QTextEdit, QPlainTextEdit, QSizePolicy, QMenu, QMessageBox
except Exception:
    # PyQt4
    from code_saturne.gui.base.QtGui import QMainWindow, QMessageBox, \
        QAction, QFileDialog, QTextEdit, QPlainTextEdit, QSizePolicy, QMenu, QMessagBox

from code_saturne.model.Common import GuiParam
#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("QCouplingEditorView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# CouplingEditorView class
#-------------------------------------------------------------------------------

_ordered_keys = ['domain', 'solver', 'n_procs_min', 'n_procs_max', 'n_procs_weight']

#-------------------------------------------------------------------------------
# domain_entry class
#-------------------------------------------------------------------------------

class domain_entry(object):
    """
    Class for domain data entry management.
    """

    # ---------------------------------------------------------------
    def __init__(self, data_dict):
        """
        Constructor
        """

        self.name   = None
        self.n_ent  = 0
        self.keys   = []
        self.labels = []
        self.vals   = []

        keys = data_dict.keys()

        for k in _ordered_keys:
            if k in keys:
                self.keys.append(k)
                self.labels.append(self.__create_label__(k))
                self.vals.append(self.__create_line__(data_dict[k]))

                if k in ('domain', 'solver'):
                    self.vals[-1].setEnabled(False)
                    if k == 'domain':
                        self.name = data_dict[k]

                self.n_ent += 1

        for k in keys:
            if k not in _ordered_keys:
                self.keys.append(k)
                self.labels.append(self.__create_label__(k))
                self.vals.append(self.__create_line__(data_dict[k]))
                self.n_ent += 1
    # ---------------------------------------------------------------

    # ---------------------------------------------------------------
    def __create_label__(self, text):

        label = QLabel()

        label.setText(text)
        label.setAlignment(Qt.AlignLeft|Qt.AlignVCenter)

        return label
    # ---------------------------------------------------------------

    # ---------------------------------------------------------------
    def __create_line__(self, val):

        line = QLineEdit()
        line.setText(str(val))

        return line
    # ---------------------------------------------------------------

#-------------------------------------------------------------------------------
# CouplingEditorView class
#-------------------------------------------------------------------------------

class CouplingEditorView(QWidget):

    # ---------------------------------------------------------------
    def __init__(self, parent, cfgfile=None):
        """
        Constructor
        """

        QWidget.__init__(self, parent)

        self.cfgfile=cfgfile

        self.data  = []
        self.pages = []

        self.data_modified = False
        self.file_loaded   = False

        self.initUI(self.cfgfile)
    # ---------------------------------------------------------------

    # ---------------------------------------------------------------
    def load_cfg_file(self, cfgfile=None):
        """
        Load a configuration file
        """

        self.cfgfile = cfgfile
        if self.cfgfile:
            from code_saturne.base import cs_run_conf
            self.run_conf = cs_run_conf.run_conf(self.cfgfile)
            if self.run_conf:
                for dom in self.run_conf.get_coupling_parameters():
                    self.add_tab_data(dom)

                if self.data != []:
                    self.file_loaded   = True
                self.data_modified = False
    # ---------------------------------------------------------------

    # ---------------------------------------------------------------
    def save_cfg_file(self):
        """
        Save the configuration file
        """

        for d in self.data:
            secname = d.name.lower()
            for entry in range(d.n_ent):
                self.run_conf.set(secname,
                                  str(d.keys[entry]),
                                  str(d.vals[entry].text()))

        self.run_conf.save()
    # ---------------------------------------------------------------

    # ---------------------------------------------------------------
    def initUI(self, cfgfile):
        """
        Create the visual layout
        """

        # Generate dictionary here
        self.load_cfg_file(cfgfile)

        if not self.file_loaded:
           title = self.tr("Warning")
           msg   = self.tr("File %s does not contain coupling parameters." % cfgfile)
           QMessageBox.warning(self, title, msg)
           return

        self.tabwidget = QTabWidget(self)
        vbox = QVBoxLayout()

        # Title
        l1 = QLabel()
        l1.setText("Coupled cases")
        l1.setAlignment(Qt.AlignCenter)
        lfont = QFont()
        lfont.setBold(True)
        l1.setFont(lfont)

        # File
        fileLabel = QLabel()
        fileLabel.setText("Parameters file: ")
        fileLabel.setFont(lfont)
        filePLabel = QLabel()
        filePLabel.setText(cfgfile)
        filePLabel.setFont(lfont)
        fileHbox = QHBoxLayout()
        fileHbox.addWidget(fileLabel)
        fileHbox.addWidget(filePLabel)
        fileHbox.addStretch(1)

        # Final layout
        vbox.addWidget(l1)
        vbox.addLayout(fileHbox)
        vbox.addWidget(self.tabwidget)

        self.setLayout(vbox)

        for i, d in enumerate(self.data):
            dom = d.name
            if dom:
                tab_name = "domain: %s" % dom
            else:
                tab_name = "domain: unknown"

            self.pages.append(self.add_tab(d))
            self.tabwidget.addTab(self.pages[i], tab_name)


        self.tabwidget.setCurrentIndex(0)
    # ---------------------------------------------------------------

    # ---------------------------------------------------------------
    def add_tab_data(self, data):
        """
        Add tab data
        """

        self.data.append(domain_entry(data))
    # ---------------------------------------------------------------

    # ---------------------------------------------------------------
    def add_tab(self, dom_data):
        """
        Add a tab to the tabWidget
        """

        rows = dom_data.n_ent
        cols = 2

        Labels = dom_data.labels
        LineEdits = dom_data.vals

        layout = QGridLayout()

        for row in range(rows):
            layout.addWidget(Labels[row], row, 0)
            layout.addWidget(LineEdits[row]  , row, 1)
            LineEdits[row].textChanged[str].connect(self.dataModifed)

        lSpacer = QLabel()
        lSpacer.setText(30*" ")
        vSpacer = QSpacerItem(40, 40, QSizePolicy.Minimum, QSizePolicy.Expanding)
        layout.addWidget(lSpacer, rows, 0)
        layout.addItem(vSpacer, rows, 1)


        tab = QWidget()
        tab.setLayout(layout)

        return tab
    # ---------------------------------------------------------------

    # ---------------------------------------------------------------
    @pyqtSlot()
    def dataModifed(self):

        self.data_modified = True
    # ---------------------------------------------------------------

#-------------------------------------------------------------------------------
# QCouplingEditor class (Independent window)
#-------------------------------------------------------------------------------

class QCouplingEditor(QMainWindow):
    """
    Independant window
    """

    # ---------------------------------------------------------------
    def __init__(self, parent=None, cfgfile=None, standalone_mode=False):
        """
        Constructor
        """

        super(QCouplingEditor, self).__init__(parent)
        self.setGeometry(50, 50, 500, 300)

        self.setWindowTitle("Coupling parameters editor")
        self.parent = parent
        self.stdalone = standalone_mode

        self.cfgfile = cfgfile

        self.editorView = CouplingEditorView(None, self.cfgfile)

        # If standalone mode, exit if file not loaded
        if self.stdalone:
            if not self.editorView.file_loaded:
                self.close()
                sys.exit(0)
                return None

        # Save file action
        save_img_path = ":/icons/22x22/document-save.png"
        icon_save     = QIcon()
        icon_save.addPixmap(QPixmap(_fromUtf8(save_img_path)),
                            QIcon.Normal,
                            QIcon.Off)
        self.saveFileAction = QAction(icon_save, "Save", self)
        self.saveFileAction.setShortcut("Ctrl+S")
        self.saveFileAction.setStatusTip('Save file')
        self.saveFileAction.triggered.connect(self.saveFile)

        # Close file action
        close_img_path = ":/icons/22x22/process-stop.png"
        icon_close     = QIcon()
        icon_close.addPixmap(QPixmap(_fromUtf8(close_img_path)),
                             QIcon.Normal,
                             QIcon.Off)
        self.closeFileAction = QAction(icon_close, "Close file", self)
        self.closeFileAction.setShortcut("Ctrl+Q")
        self.closeFileAction.setStatusTip('Close opened file')
        self.closeFileAction.triggered.connect(self.closeOpenedFile)

        # Toolbar
        self.toolbar = self.addToolBar("Options")
        self.toolbar.addAction(self.saveFileAction)
        self.toolbar.addAction(self.closeFileAction)

        self.setCentralWidget(self.editorView)
    # ---------------------------------------------------------------

    # ---------------------------------------------------------------
    def saveFile(self):
        """
        Save the configuration file
        """

        self.editorView.save_cfg_file()
        self.editorView.data_modified = False
    # ---------------------------------------------------------------

    # ---------------------------------------------------------------
    def closeOpenedFile(self):
        """
        Close and save if necessary the opened configuration file
        """

        if self.editorView.file_loaded:
            choice = QMessageBox.question(self, 'Coupling parameters editor',
                                          'Exit editor ?',
                                          QMessageBox.Yes | QMessageBox.No)
        else:
            choice = QMessageBox.Yes

        if choice == QMessageBox.Yes:

            if self.editorView.data_modified:
                choice = QMessageBox.question(self, 'Coupling parameters editor',
                                              'Save modifications ?',
                                              QMessageBox.Yes | QMessageBox.No)
                if choice == QMessageBox.Yes:
                    self.saveFile()
                else:
                    pass

            settings = QSettings()
            settings.setValue("MainWindow/Geometry",
                              self.saveGeometry())

            self.editorView.file_loaded = False

            self.close()
            return 0
        else:
            return 1
    # ---------------------------------------------------------------

#-------------------------------------------------------------------------------
# END OF FILE
#-------------------------------------------------------------------------------
