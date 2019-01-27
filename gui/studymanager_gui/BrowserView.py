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
This module defines the following classes:
- BrowserView
"""

#-------------------------------------------------------------------------------
# Standard modules
#-------------------------------------------------------------------------------

import sys, logging


#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------

from code_saturne.Base.QtCore    import *
from code_saturne.Base.QtGui     import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.Base.BrowserForm import Ui_BrowserForm
from code_saturne.model.Common import GuiParam
from code_saturne.Base.BrowserView import TreeItem, TreeModel
from code_saturne.Base.BrowserView import BrowserView as SaturneBrowserView
from code_saturne.studymanager_gui.Toolbox import displaySelectedPage

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("BrowserView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

class BrowserView(SaturneBrowserView):
    """
    Class for the browser widget
    """
    def __init__(self):
        """
        Constructor
        """
        SaturneBrowserView.__init__(self)
        tree = self._browser()
        self.model = TreeModel(str(tree))

        self.treeView.setModel(self.model)

    def _browser(self):
        tree ="""
Paths
Manage cases
Define plotter
"""
        return tree


    def activeSelectedPage(self, index):
        """
        """
        self.treeView.selectionModel().select(index, QItemSelectionModel.SelectCurrent)

        return


    def display(self, root, case, stbar, tree):
        """
        """
        index = self.treeView.currentIndex()
        item  = index.internalPointer()
        name  = item.itemData[0]
        case['current_tab'] = 0
        case['current_index'] = index
        return displaySelectedPage(name, root, case, stbar, tree)


    def configureTree(self, case):
        """
        Public method.
        Configures the browser with users data.
        """
        self.case = case

        self.setRowOpen(self.tr('Paths'))
        self.setRowOpen(self.tr('Manage cases'))
        self.setRowOpen(self.tr('Define plotter'))

        self.__hideRow()


#-------------------------------------------------------------------------------
# Testing part
#-------------------------------------------------------------------------------

if __name__ == "__main__":
    app = QApplication(sys.argv)
    BrowserView = BrowserView()
    BrowserView.show()
    sys.exit(app.exec_())

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
