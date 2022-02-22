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

from code_saturne.gui.base.QtCore    import *
from code_saturne.gui.base.QtGui     import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.gui.base.BrowserForm import Ui_BrowserForm
from code_saturne.model.Common import GuiParam
from code_saturne.gui.base.BrowserView import TreeItem, TreeModel
from code_saturne.gui.base.BrowserView import BrowserView as SaturneBrowserView
from code_saturne.gui.studymanager_gui.Toolbox import displaySelectedPage

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
    def __init__(self, parent=None):
        """
        Constructor
        """
        self.parent = parent

        SaturneBrowserView.__init__(self, parent)
        tree = self._browser()
        self.model = TreeModel(tree)

        self.treeView.setModel(self.model)

    def _browser(self):
        tree = [TreeItem(['Pages'], 'folder')]

        for section in ("Paths", "Manage cases", "Define plotter"):
            section_item = TreeItem([section], typename='folder-new', parent=tree[0])
            tree[0].appendChild(section_item)

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
