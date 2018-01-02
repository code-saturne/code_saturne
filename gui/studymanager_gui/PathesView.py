# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2018 EDF S.A.
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
This module defines the Page in which the user defines the path of the treated
case. The Page verify that the case configuration directories is appropriate.

This module contains the following classes:
- PathesView
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import os, string
import logging

from code_saturne.Base.QtCore    import *
from code_saturne.Base.QtGui     import *
from code_saturne.Base.QtWidgets import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.Base.Toolbox import GuiParam
from code_saturne.studymanager_gui.PathesForm import Ui_PathesForm
from code_saturne.studymanager_gui.PathesModel import PathesModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("PathesView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

class PathesView(QWidget, Ui_PathesForm):
    """
    Class to open Pathes Page.
    """
    def __init__(self, parent, case):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_PathesForm.__init__(self)
        self.setupUi(self)

        self.case = case
        self.parent = parent

        self.mdl = PathesModel(self.case)

        # Create the Page layout.
        self.toolButtonRepository.pressed.connect(self.searchDirRepo)
        self.toolButtonDestination.pressed.connect(self.searchDirDest)
        self.lineEditRepository.setReadOnly(True)
        self.lineEditDestination.setReadOnly(True)

        self.repo_path = self.mdl.getRepositoryPath()
        self.lineEditRepository.setText(self.repo_path)
        self.dest_path = self.mdl.getDestinationPath()
        self.lineEditDestination.setText(self.dest_path)

        if self.repo_path == "":
            self.toolButtonRepository.setStyleSheet("background-color: red")
        else:
            self.toolButtonRepository.setStyleSheet("background-color: green")
        if self.dest_path == "":
            self.toolButtonDestination.setStyleSheet("background-color: red")
        else:
            self.toolButtonDestination.setStyleSheet("background-color: green")


    def searchDirRepo(self):
        """
        Open a File Dialog in order to search the repository directory.
        """
        if self.repo_path != "":
            path_case = self.repo_path
        else:
            path_case = os.path.dirname(os.getcwd())

        title    = self.tr("Select directory")
        default  = path_case
        options  = QFileDialog.ShowDirsOnly # | QFileDialog.DontResolveSymlinks
        dir_name = QFileDialog.getExistingDirectory(self, title, default, options)
        dir_name = str(dir_name)

        if dir_name:
            self.repo_path = dir_name
            self.lineEditRepository.setText(self.repo_path)
            self.toolButtonRepository.setStyleSheet("background-color: green")
            self.mdl.setRepositoryPath(self.repo_path)

        if self.mdl.getStudyNumber() == 0:
            title = self.tr("Loading study/case")
            msg   = self.tr("Do you want to use automatic load ?")

            reply = QMessageBox.question(self,
                                         title,
                                         msg,
                                         QMessageBox.Yes|
                                         QMessageBox.No)
            if reply == QMessageBox.Yes:
                self.mdl.loadStudyAndCases(self.repo_path)


    def searchDirDest(self):
        """
        Open a File Dialog in order to search the destination directory.
        """
        if self.dest_path != "":
            path_case = self.dest_path
        else:
            path_case = os.path.dirname(os.getcwd())

        title    = self.tr("Select directory")
        default  = path_case
        options  = QFileDialog.ShowDirsOnly # | QFileDialog.DontResolveSymlinks
        dir_name = QFileDialog.getExistingDirectory(self, title, default, options)
        dir_name = str(dir_name)

        if dir_name:
            self.dest_path = dir_name
            self.lineEditDestination.setText(self.dest_path)
            self.toolButtonDestination.setStyleSheet("background-color: green")
            self.mdl.setDestinationPath(self.dest_path)


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
