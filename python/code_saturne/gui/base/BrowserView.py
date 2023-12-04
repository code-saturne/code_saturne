# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2023 EDF S.A.
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
#Third-party modules
#-------------------------------------------------------------------------------

from code_saturne.gui.base.QtCore import *
from code_saturne.gui.base.QtGui import *
from code_saturne.gui.base.QtWidgets import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.model.LocalizationModel import LocalizationModel
from code_saturne.model.Common import GuiParam

from code_saturne.gui.base.BrowserForm import Ui_BrowserForm
from code_saturne.gui.base.Toolbox import displaySelectedPage
from code_saturne.gui.base.QtPage import from_qvariant, to_text_string

import code_saturne.gui.base.resource_base_rc

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("BrowserView")
log.setLevel(GuiParam.DEBUG)

try:
    _fromUtf8 = QString.fromUtf8
except Exception:
    def _fromUtf8(s):
        return s


class TreeItem:
    """ Reimplementation of the Qt TreeItem.
        Most attribute and method names are not PEP-8 compliant but follow the Qt convention.
        See : https://doc.qt.io/qt-5/qtwidgets-itemviews-simpletreemodel-example.html
        """

    def __init__(self, data, typename, parent=None):
        self.parentItem = parent
        self.itemData = data
        self.itemType = typename
        self.itemIcon = None
        self.childItems = []
        if parent is None:
            self.n_parents = 0
        else:
            self.n_parents = parent.n_parents + 1

    def appendChild(self, item):
        self.childItems.append(item)

    def insertChildren(self, position, count, columns):
        for row in range(count):
            data = [None for v in range(columns)]
            item = TreeItem(data, typename="folder-new", parent=self)
            self.childItems.insert(position, item)
        return True

    def removeChildren(self, position, count):
        if (position < 0) or (position + count > len(self.childItems)):
            return False
        for row in range(count):
            self.childItems.pop(position)

    def child(self, row):
        return self.childItems[row]

    def childWithName(self, name):
        for child in self.childItems:
            if child.data(0) == name:
                return child
        return None

    def remove(self):
        while self.childCount() > 0:
            self.child(-1).remove()

        self.parentItem.childItems.remove(self)
        return

    def childCount(self):
        return len(self.childItems)

    def columnCount(self):
        return len(self.itemData)

    def data(self, column):
        return self.itemData[column]

    def setData(self, column, value):
        self.itemData[column] = value
        return True

    def parent(self):
        return self.parentItem

    def row(self):
        if self.parentItem:
            return self.parentItem.childItems.index(self)

        return 0

    def level(self):
        return self.n_parents


#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

class TreeModel(QAbstractItemModel):
    """A model representing the widget tree structure.
    """

    def __init__(self, data, parent=None):
        """Constructs a new item model with the given I{parent}.

        @type data: C{QString}
        @param data: content of the new item
        @type parent: C{QObject} or C{None}
        @param parent: parent of the new item
        """
        QAbstractItemModel.__init__(self, parent)

        self.rootItem = data[0]

    def columnCount(self, parent):
        """Returns the number of columns for the children of the given I{parent}.

        @type parent: C{QModelIndex}
        @param parent: parent of the item
        @return: C{int}
        """
        if parent.isValid():
            return parent.internalPointer().columnCount()
        else:
            return self.rootItem.columnCount()

    def data(self, index, role):
        """Returns the data stored under the given I{role} for the item referred to by the I{index}.

        @type index: C{QModelIndex}
        @param index: used to locate data in a data model
        @type role: C{Qt.ItemDataRole}
        @param role: used by the view to indicate to the model which type of data it needs
        @return: C{QVariant}
        """
        if not index.isValid():
            return None

        item = index.internalPointer()
        column = index.column()

        if role == Qt.DisplayRole:
            # return text for columns
            if column == 0:
                return item.itemData[column]

        elif role == Qt.DecorationRole:
            # return icon for first column
            if column == 0:
                page_name = item.itemData[0]

                style = QWidget().style()
                icon = QIcon()

                if page_name == self.tr('Mesh'):
                    img_path = ":/icons/22x22/cube_mesh.png"
                    icon.addPixmap(QPixmap(_fromUtf8(img_path)), QIcon.Normal, QIcon.Off)
                elif page_name == self.tr('Calculation features'):
                    img_path = ":/icons/22x22/calculation_features.png"
                    icon.addPixmap(QPixmap(_fromUtf8(img_path)), QIcon.Normal, QIcon.Off)
#                elif page_name == self.tr('Volume zones'):
#                    img_path = ":/icons/22x22/volume_zones.png"
#                    icon.addPixmap(QPixmap(_fromUtf8(img_path)), QIcon.Normal, QIcon.Off)
                elif page_name == self.tr('Volume conditions'):
                    img_path = ":/icons/22x22/cube_volume_zone_v2.png"
                    icon.addPixmap(QPixmap(_fromUtf8(img_path)), QIcon.Normal, QIcon.Off)
#                elif page_name == self.tr('Boundary zones'):
#                    img_path = ":/icons/22x22/boundary_conditions.png"
#                    icon.addPixmap(QPixmap(_fromUtf8(img_path)), QIcon.Normal, QIcon.Off)
                elif page_name == self.tr('Boundary conditions'):
                    img_path = ":/icons/22x22/cube_bc.png"
                    icon.addPixmap(QPixmap(_fromUtf8(img_path)), QIcon.Normal, QIcon.Off)
                elif page_name == self.tr('Time settings'):
                    img_path = ":/icons/22x22/time_stepping.png"
                    icon.addPixmap(QPixmap(_fromUtf8(img_path)), QIcon.Normal, QIcon.Off)
                elif page_name == self.tr('Numerical parameters'):
                    img_path = ":/icons/22x22/numerical_params.png"
                    icon.addPixmap(QPixmap(_fromUtf8(img_path)), QIcon.Normal, QIcon.Off)
                elif page_name == self.tr('Postprocessing'):
                    img_path = ":/icons/22x22/postprocessing.png"
                    icon.addPixmap(QPixmap(_fromUtf8(img_path)), QIcon.Normal, QIcon.Off)
                elif page_name == 'Closure modeling':
                    img_path = ":/icons/22x22/closure_modeling.png"
                    icon.addPixmap(QPixmap(_fromUtf8(img_path)), QIcon.Normal, QIcon.Off)
                elif page_name == self.tr('Performance settings'):
                    img_path = ":/icons/22x22/run_parameters.png"
                    icon.addPixmap(QPixmap(_fromUtf8(img_path)), QIcon.Normal, QIcon.Off)

                elif item.itemType == "folder-new":
                    icon = style.standardIcon(QStyle.SP_FileLinkIcon)
                elif item.itemType == "folder-close":
                    icon = style.standardIcon(QStyle.SP_FileIcon)
                elif item.itemType == "folder-open":
                    icon = style.standardIcon(QStyle.SP_FileIcon)
                elif item.itemType == "file-open":
                    icon = style.standardIcon(QStyle.SP_FileIcon)
                elif item.itemType == "file-new":
                    icon = style.standardIcon(QStyle.SP_FileLinkIcon)
                if sys.platform.startswith("win"):
                    if item.itemType == "file-open":
                        icon = style.standardIcon(QStyle.SP_ToolBarHorizontalExtensionButton)
                    elif item.itemType == "file-new":
                        icon = style.standardIcon(QStyle.SP_ToolBarHorizontalExtensionButton)

                return icon

        return None

    def flags(self, index):
        """What we can do with the item.

        @type index: C{QModelIndex}
        @param index: used to locate data in a data model
        """
        if not index.isValid():
            return Qt.ItemIsEnabled

        flags = Qt.ItemIsEnabled | Qt.ItemIsSelectable

        return flags

    def headerData(self, section, orientation, role):
        """Return the header of the tree.*

        @return: C{QVariant}
        """
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return self.rootItem.data(section)

        return None

    def index(self, row, column, parent):
        """Returns the index of the item in the model specified by the given I{row}, I{column} and I{parent} index.

        @type row: C{int}
        @param row: row of the item
        @type column: C{int}
        @param column: column of the item
        @type parent: C{QModelIndex}
        @param parent: parent of the item
        @return: C{QModelIndex}
        """
        if not parent.isValid():
            parentItem = self.rootItem
        else:
            parentItem = parent.internalPointer()

        # FIXME: why childItem can be None?
        try:
            childItem = parentItem.child(row)
        except:
            childItem = None

        if childItem:
            return self.createIndex(row, column, childItem)
        else:
            return QModelIndex()


    def parent(self, index):
        """Returns the parent of the model item with the given index.

        @type index: C{QModelIndex}
        @param index: index of the child
        @return: C{QModelIndex}
        """
        if not index.isValid():
            return QModelIndex()

        childItem = index.internalPointer()
        parentItem = childItem.parent()

        if parentItem == self.rootItem:
            return QModelIndex()

        return self.createIndex(parentItem.row(), 0, parentItem)

    def rowCount(self, parent):
        """Returns the number of rows under the given I{parent}.

        @type parent: C{QModelIndex}
        @param parent: parent of the item
        @return: C{int}
        """

        if not parent.isValid():
            parentItem = self.rootItem
        else:
            parentItem = parent.internalPointer()

        return parentItem.childCount()

    def match(self, start, role, value, hits, flags):
        """
        @type start: C{QModelIndex}
        @type role: C{Qt.ItemDataRole}
        @type value: C{QVarient}
        @type hits: C{int}
        @type flags: C{Qt.MatchFlag}
        """
        result = []
        p = self.parent(start)
        st = start.row()
        to = self.rowCount(p)

        for r in range(st, to):

            index = self.index(r, start.column(), p)
            if not index.isValid():
                pass
            v = self.data(index, role)

            if flags == Qt.MatchExactly:
                if value == v:
                    result.append(index)
            else:
                continue

            if self.hasChildren(index):
                result += self.match(self.index(0, index.column(), index),
                                     role, value, hits, flags)

        return result

    def itemLocalization(self, data, role=Qt.DisplayRole):
        info = []
        search_item = from_qvariant(data, to_text_string)

        start = self.index(0, 0, QModelIndex())
        indexList = self.match(start, role, search_item, -1, Qt.MatchExactly)

        for index in indexList:
            item = index.internalPointer()
            column = index.column()
            row = index.row()
            parent = self.parent(index)

            info.append((row, column, parent))

        return info

    def appendRowWithValue(self, value, parent=QModelIndex()):
        nb_items = len(self.getItem(parent).childItems)
        self.insertRows(nb_items, 1, parent)
        last_index = self.index(nb_items, 0, parent)
        self.setData(last_index, value, role=Qt.EditRole)
        return

    def insertRows(self, position, rows, parent=QModelIndex()):
        parentItem = self.getItem(parent)
        self.beginInsertRows(parent, position, position + rows - 1)
        success = parentItem.insertChildren(position, rows, self.rootItem.columnCount())
        self.endInsertRows()
        return success

    def removeRows(self, position, rows, parent=QModelIndex()):
        parentItem = self.getItem(parent)
        self.beginRemoveRows(parent, position, position + rows - 1)
        success = parentItem.removeChildren(position, rows)
        self.endRemoveRows()
        return success

    def setData(self, index, value, role=Qt.EditRole):
        if role != Qt.EditRole:
            return False
        item = self.getItem(index)
        result = item.setData(index.column(), value)
        if result:
            self.dataChanged.emit(index, index)
        return result

    def getItem(self, index):
        if index.isValid():
            item = index.internalPointer()
            if item:
                return item
        return self.rootItem


class BrowserView(QWidget, Ui_BrowserForm):
    """
    Class for the browser widget
    """

    def __init__(self, parent=None):
        """
        Constructor
        """
        self.parent = parent

        QWidget.__init__(self, parent)

        Ui_BrowserForm.__init__(self)
        self.setupUi(self)

        tree = self._browser()
        self.model = TreeModel(tree)
        # TODO check how self.model should be initialized
        # self.model = TreeModel(from_qvariant(tree, to_text_string))

        self.treeView.setModel(self.model)
        self.treeView.header().hide()
        self.treeView.setAnimated(True)
        self.treeView.setVerticalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        self.treeView.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
        self.treeView.setAlternatingRowColors(True)
        self.treeView.setWindowTitle("Simple Tree Model")

        # Popup menu
        self.treeView.setContextMenuPolicy(Qt.CustomContextMenu)
        self.treeView.customContextMenuRequested[QPoint].connect(self.displayPopup)

        # Receive change in selection
        self.treeView.pressed[QModelIndex].connect(self.onItemPressed)
        self.treeView.expanded[QModelIndex].connect(self.onFolderOpen)
        self.treeView.collapsed[QModelIndex].connect(self.onFolderClose)

    def _getSectionList(self):

        sl = ['Calculation environment', 'Mesh', 'Calculation features',
              'Closure modeling',
              'Particles and droplets tracking', 'Volume conditions',
              'Boundary conditions', 'Coupling parameters', 'Time settings',
              'Numerical parameters',
              'Postprocessing', 'Performance settings']

        return sl

    def _getSubsectionList(self, section):

        if section == 'Calculation environment':
            return ['Notebook', 'Time tables']
        elif section == 'Mesh':
            return ['Preprocessing', "Volume zones", "Boundary zones"]
        elif section == 'Calculation features':
            return ['Main fields', 'Deformable mesh', 'Gas combustion',
                    'Pulverized fuel combustion', 'Electrical models',
                    'Atmospheric flows', 'Thermal model', 'Body forces',
                    'Turbulence models', 'Species transport',
                    'Groundwater flows', 'Turbomachinery',
                    'Fans', 'Non condensable gases']
        elif section == 'Closure modeling':
            return ['Interfacial area',
                    'Interfacial enthalpy transfer',
                    'Wall transfer parameters',
                    'Particles interactions']
        elif section == 'Particles and droplets tracking':
            return ['Statistics']
        elif section == 'Volume conditions':
            return []
        elif section == "Boundary conditions":
            return []
        elif section == 'Coupling parameters':
            return ['Immersed Boundaries']
        elif section == 'Time settings':
            return ['Start/Restart']
        elif section == 'Numerical parameters':
            return ['Equation parameters']
        elif section == 'Postprocessing':
            return ['Calculator', 'Additional user arrays', 'Time averages',
                    'Volume solution control', 'Surface solution control',
                    'Lagrangian solution control', 'Profiles',
                    'Balance by zone']
        elif section == 'Performance settings':
            return []
        else:
            return None

    def _browser(self):

        # TODO : why is tree a list ?
        tree = [TreeItem(['Pages'], 'folder')]

        for section in self._getSectionList():
            section_item = TreeItem([section], typename='folder-new', parent=tree[0])
            tree[0].appendChild(section_item)

            for subsection in self._getSubsectionList(section):
                subsection_item = TreeItem([subsection], typename='file-new', parent=section_item)
                section_item.appendChild(subsection_item)

        return tree

    def setRowClose(self, string):
        log.debug("setRowClose(): %s" % string)
        item_info_list = self.model.itemLocalization(string)
        for item_info in item_info_list:
            row = item_info[0]
            parent = item_info[2]
            self.treeView.setRowHidden(row, parent, True)

    def setRowOpen(self, string):
        log.debug("setRowOpen(): %s" % string)
        item_info_list = self.model.itemLocalization(string)
        for item_info in item_info_list:
            row = item_info[0]
            parent = item_info[2]
            self.treeView.setRowHidden(row, parent, False)

    def setRowShow(self, string, status=True):
        log.debug("setRowVisible(): %s" % string)
        item_info_list = self.model.itemLocalization(string)
        hidden = not status
        for item_info in item_info_list:
            row = item_info[0]
            parent = item_info[2]
            self.treeView.setRowHidden(row, parent, hidden)

    def isRowClose(self, string):
        log.debug("isRowClose(): %s" % string)
        item_info_list = self.model.itemLocalization(string)
        for item_info in item_info_list:
            row = item_info[0]
            column = item_info[1]
            parent = item_info[2]
            index = self.model.index(row, column, parent)
            # FIXME: this return should not be in a loop
            return self.treeView.isRowHidden(row, index)

    @pyqtSlot('QModelIndex')
    def onItemPressed(self, index):
        item = index.internalPointer()
        if item.itemType == "file-new":
            item.itemType = "file-open"

    @pyqtSlot('QModelIndex')
    def onFolderOpen(self, index):
        """
        public slot

        change the item type when the folder is opened

        @type index: C{QModelIndex}
        @param index: index in the model of the selected folder
        """
        item = index.internalPointer()
        if item.itemType == "folder-new" or item.itemType == "folder-close":
            item.itemType = "folder-open"

    @pyqtSlot('QModelIndex')
    def onFolderClose(self, index):
        """
        public slot

        change the item type when the folder is closed

        @type index: C{QModelIndex}
        @param index: index in the model of the selected folder
        """
        item = index.internalPointer()
        if item.itemType == "folder-new" or item.itemType == "folder-open":
            item.itemType = "folder-close"

    @pyqtSlot()
    def displayPopup(self):
        """
        public slot

        create the popup menu of the Browser
        """
        self.fileMenu = QMenu(self.treeView)

        self.actionExpand = QAction(self.tr("Expand"), self.treeView)
        self.actionExpand.triggered.connect(self.openTreeFolder)

        self.actionCollapse = QAction(self.tr("Collapse"), self.treeView)
        self.actionCollapse.triggered.connect(self.closeTreeFolder)

        # ... TODO
        # self.actionWelcome = QAction(self.tr("Welcome page"), self.treeView)

        self.fileMenu.addAction(self.actionExpand)
        self.fileMenu.addAction(self.actionCollapse)

        cursor = QCursor()
        self.fileMenu.popup(cursor.pos())
        self.fileMenu.show()

    def activeSelectedPage(self, index):
        """
        """
        if index != None:
            self.treeView.selectionModel().select(index,
                                                  QItemSelectionModel.SelectCurrent)

        return

    def addBrowserZone(self, zone_name, parent_name):
        itemInfoList = self.model.itemLocalization(parent_name)
        for itemInfo in itemInfoList:
            parent = self.model.index(*itemInfo)
            self.model.appendRowWithValue(zone_name, parent)

    def removeAllBrowserZones(self, parent_name):
        localization = self.model.itemLocalization(parent_name)
        parent_index = self.model.index(*localization[0])
        parent_item = self.model.getItem(parent_index)
        nb_zones = parent_item.childCount()
        self.model.removeRows(0, nb_zones, parent_index)

    def updateBrowserZones(self, new_zone_names, parent_name):
        self.removeAllBrowserZones(parent_name)
        for zone_name in new_zone_names:
            self.addBrowserZone(zone_name, parent_name)

    def display(self, root, case, stbar, tree):
        """
        """
        index = self.treeView.currentIndex()
        item = index.internalPointer()
        name = item.itemData[0]
        case['current_tab'] = 0
        case['current_index'] = index
        return displaySelectedPage(name, root, case, stbar, tree)

    def isFolder(self):
        """
        Return True if current item is a folder (parent)
        """
        index = self.treeView.currentIndex()
        item = index.internalPointer()
        return item.childCount() != 0

    def openSingleFolder(self, string):
        """
        Open a single folder of the Tree.
        """
        itemInfoList = self.model.itemLocalization(string)
        for itemInfo in itemInfoList:
            row = itemInfo[0]
            column = itemInfo[1]
            parent = itemInfo[2]
            index = self.model.index(row, column, parent)
            self.treeView.expand(index)

    @pyqtSlot()
    def openTreeFolder(self):
        """
        public slot

        open all folders of the Tree.
        """
        self.treeView.expandAll()

        parent = QModelIndex()
        column = 0
        for row in range(self.model.rowCount(parent)):
            index = self.model.index(row, column, parent)
            self.onFolderOpen(index)

        if hasattr(self, 'case'):
            self.configureTree(self.case)

    def closeSingleFolder(self, string):
        """
        Close a single folder of the Tree.
        """
        itemInfoList = self.model.itemLocalization(string)
        for itemInfo in itemInfoList:
            row = itemInfo[0]
            column = itemInfo[1]
            parent = itemInfo[2]
            index = self.model.index(row, column, parent)
            self.treeView.collapse(index)

    @pyqtSlot()
    def closeTreeFolder(self):
        """
        public slot

        close all folders of the Tree.
        """
        self.treeView.collapseAll()

        parent = QModelIndex()
        column = 0
        for row in range(self.model.rowCount(parent)):
            index = self.model.index(row, column, parent)
            self.onFolderClose(index)

    def __configureTreePrepro(self, case):
        """
        Public method.
        Configures the browser with users data.
        """

        show_prepro = True
        if case['run_type'] == 'none':
            show_prepro = False

        self.setRowShow(self.tr('Preprocessing'), show_prepro)
        self.setRowShow(self.tr('Postprocessing'), show_prepro)

        self.setRowClose(self.tr('Calculation features'))
        """
        self.setRowClose(self.tr('Main fields'))
        self.setRowClose(self.tr('Deformable mesh'))
        self.setRowClose(self.tr('Gas combustion'))
        self.setRowClose(self.tr('Pulverized fuel combustion'))
        self.setRowClose(self.tr('Electrical models'))
        self.setRowClose(self.tr('Atmospheric flows'))
        self.setRowClose(self.tr('Thermal model'))
        self.setRowClose(self.tr('Turbulence models'))
        self.setRowClose(self.tr('Transported Variables'))
        self.setRowClose(self.tr('Turbomachinery'))
        self.setRowClose(self.tr('Groundwater flows'))
        self.setRowClose(self.tr('Fans'))
        self.setRowClose(self.tr('Non condensable gases'))
        """

        self.setRowClose(self.tr('Closure modeling'))
        """
        self.setRowClose(self.tr('Interfacial area'))
        self.setRowClose(self.tr('Interfacial enthalpy transfer'))
        self.setRowClose(self.tr('Wall transfer parameters'))
        self.setRowClose(self.tr('Particles interactions'))
        """

#        """
#        self.setRowClose(self.tr('Body forces'))
#        """

        self.setRowClose(self.tr('Particles and droplets tracking'))
        """
        self.setRowClose(self.tr('Statistics'))
        """

        self.setRowShow(self.tr('Volume conditions'), False)
        """
        self.setRowClose(self.tr('Initialization'))
        self.setRowClose(self.tr('Head losses'))
        self.setRowClose(self.tr('Porosity'))
        self.setRowClose(self.tr('Source terms'))
        self.setRowClose(self.tr('Groundwater laws'))
        """
        """
        self.setRowClose(self.tr('Boundary_conditions'))
        self.setRowClose(self.tr('Immersed Boundaries'))
        """

        self.setRowClose(self.tr('Time settings'))

        self.setRowClose(self.tr('Numerical parameters'))
        """
        self.setRowClose(self.tr('Equation parameters'))
        """

        self.setRowClose(self.tr('Time averages'))
        self.setRowClose(self.tr('Calculator'))
        self.setRowClose(self.tr('Additional user arrays'))
        self.setRowClose(self.tr('Volume solution control'))
        self.setRowClose(self.tr('Surface solution control'))
        self.setRowClose(self.tr('Lagrangian solution control'))
        self.setRowClose(self.tr('Profiles'))
        self.setRowClose(self.tr('Balance by zone'))

        self.__hideRow()

    def configureTree(self, case):
        """
        Public method.
        Configures the browser with users data.
        """

        if case['run_type'] != 'standard':
            return self.__configureTreePrepro(case)

        p_module = ''
        if case.xmlRootNode().tagName == "NEPTUNE_CFD_GUI":
            p_module = 'neptune_cfd'

        # Precompute some values
        # ----------------------

        m_tbm = False
        m_fans = False
        m_lagr = False
        m_ale = False
        m_thermal = 0
        m_cht = False
        m_gas_comb = False
        m_sf_comb = False
        m_elec = False
        m_atmo = False
        m_gwf = False
        m_icm = False

        node_pm = case.xmlGetNode('thermophysical_models')

        if node_pm:

            node = node_pm.xmlGetNode('ale_method')
            if node and node['status'] == 'on':
                m_ale = True
            node = node_pm.xmlGetNode('turbomachinery', 'model')
            if node and node['model'] != "off":
                m_tbm = True
            node = node_pm.xmlGetNode('fans')
            if node:
                m_fans = True

            if not m_thermal:
                node = node_pm.xmlGetNode('gas_combustion', 'model')
                if node and node['model'] in ('ebu', 'd3p', 'lwp'):
                    m_gas_comb = True
                    m_thermal = 1
            if not m_thermal:
                node = node_pm.xmlGetNode('solid_fuels', 'model')
                if node and node['model'] in ('homogeneous_fuel',
                                              'homogeneous_fuel_moisture',
                                              'homogeneous_fuel_moisture_lagr'):
                    m_sf_comb = True
                    m_thermal = 1
            if not m_thermal:
                node = node_pm.xmlGetNode('joule_effect', 'model')
                if node and node['model'] in ('joule', 'arc'):
                    m_elec = True
                    m_thermal = 1
            if not m_thermal:
                node = node_pm.xmlGetNode('atmospheric_flows', 'model')
                if node and node['model'] != 'off':
                    m_atmo = True
                    if node['model'] == 'constant':
                        m_thermal = -1
                    elif node['model'] in ('dry', 'humid', 'humid_ctwr'):
                        m_thermal = 1
            if not m_thermal:
                node = node_pm.xmlGetNode('compressible_model', 'model')
                if node and node['model'] != 'off':
                    m_cpr = True
                    m_thermal = 1

            if not m_thermal:
                node = node_pm.xmlGetNode('thermal_scalar', 'model')
                if node and node['model'] != 'off':
                    m_thermal = 2

            node = node_pm.xmlGetNode('groundwater_model', 'model')
            if node and node['model'] != 'off':
                m_gwf = True
                m_thermal = -1
                m_fans = False

            if m_thermal > 0:
                m_cht = True

            if not m_icm:
                node = node_pm.xmlGetChildNode("internal_coupling")
                if node:
                    node = node.xmlGetChildNode("solid_zones")
                    if node:
                        if len(node.xmlGetChildNodeList("zone")) > 0:
                            m_icm = True
                        else:
                            m_icm = False

        node = case.xmlGetNode('lagrangian', 'model')
        if node and node['model'] != "off":
            m_lagr = True

        # Options for NEPTUNE_CFD if present:

        m_ncfd = {}
        m_ncfd['non_condens'] = False
        m_ncfd['particles_interactions'] = False
        m_ncfd['itf_area'] = False
        m_ncfd['itf_h_transfer'] = False
        m_ncfd['wall_transfer'] = False

        ncfd_fields = 0

        if p_module == 'neptune_cfd' and node_pm:
            m_fans = False
            m_thermal = 1
            m_cht = True

            node_f = node_pm.xmlInitNode('fields')
            fields = node_f.xmlGetNodeList('field')
            ncfd_fields = len(fields)

        if ncfd_fields > 1:
            from code_saturne.model.MainFieldsModel import MainFieldsModel
            predefined_flow = MainFieldsModel(case).getPredefinedFlow()
            phase_change_transfer = MainFieldsModel(case).getPhaseChangeTransferStatus()

            if (len(MainFieldsModel(case).getSolidPhaseList()) > 0):
                m_ncfd['particles_interactions'] = True

            m_ncfd['itf_area'] = True
            m_ncfd['non_condens'] = True
            m_ncfd["itf_h_transfer"] = {"on": True, "off": False}[phase_change_transfer] or m_ncfd['particles_interactions']
            m_ncfd['wall_transfer'] = {"on": True, "off": False}[phase_change_transfer]

            if predefined_flow == "free_surface":
                m_ncfd['itf_area'] = False
            elif predefined_flow == "particles_flow":
                m_ncfd['particles_interactions'] = True
                m_ncfd['wall_transfer'] = False
            elif predefined_flow == "None":
                m_ncfd['particles_interactions'] = True

        is_ncfd = (p_module == 'neptune_cfd')

        # Manage visibility
        # -----------------

        self.setRowOpen(self.tr('Notebook'))

        # Preprocessing

        self.setRowShow(self.tr('Preprocessing'), True)

        # Thermophysical Models

        self.setRowShow(self.tr('Calculation features'), True)
        self.setRowShow(self.tr('Main fields'), (p_module == 'neptune_cfd'))
        self.setRowShow(self.tr('Deformable mesh'), m_ale)
        self.setRowShow(self.tr('Gas combustion'), m_gas_comb)
        self.setRowShow(self.tr('Pulverized fuel combustion'), m_sf_comb)
        self.setRowShow(self.tr('Electrical models'), m_elec)
        self.setRowShow(self.tr('Atmospheric flows'), m_atmo)
        self.setRowShow(self.tr('Thermal model'), (m_thermal > -1))
        self.setRowShow(self.tr('Body forces'), True)
        self.setRowShow(self.tr('Turbulence models'), not m_gwf)
        self.setRowShow(self.tr('Transported Variables'))
        self.setRowShow(self.tr('Groundwater flows'), m_gwf)
        self.setRowShow(self.tr('Turbomachinery'), m_tbm)
        self.setRowShow(self.tr('Fans'), m_fans)

        self.setRowShow(self.tr('Non condensable gases'), m_ncfd['non_condens'])

        # Closure modeling

        self.setRowShow(self.tr('Closure modeling'), (ncfd_fields > 1))
        self.setRowShow(self.tr('Interfacial enthalpy transfer'), m_ncfd['itf_h_transfer'])
        self.setRowShow(self.tr('Interfacial area'), m_ncfd['itf_area'])
        self.setRowShow(self.tr('Wall transfer parameters'), m_ncfd['wall_transfer'])
        self.setRowShow(self.tr('Particles interactions'), m_ncfd['particles_interactions'])

        # Particles and droplets tracking

        self.setRowShow(self.tr('Particles and droplets tracking'), m_lagr)
        self.setRowShow(self.tr('Statistics'), m_lagr)

        # Volume zones

        self.setRowShow(self.tr('Volume zones'), True)
        self.setRowShow(self.tr('Volume conditions'), True)

        node_domain = case.xmlGetNode('solution_domain')

        # Boundary zones

        self.setRowShow(self.tr('Boundary conditions'))
        # Immersed boundaries is deactivated for the moment. Will be
        # reactivated following v6.1 once Page is updated in NCFD
        self.setRowShow(self.tr('Immersed Boundaries'), False)
        self.setRowShow(self.tr("Coupling parameters"), m_lagr or m_ale or is_ncfd or m_cht or m_icm)

        # Time settings

        self.setRowShow(self.tr('Time settings'))

        # Numerical parameters

        self.setRowShow(self.tr('Numerical parameters'))
        self.setRowShow(self.tr('Equation parameters'))

        # Postprocessing

        self.setRowShow(self.tr('Postprocessing'), True)
        self.setRowShow(self.tr('Time averages'), True)
        self.setRowShow(self.tr('Calculator'))
        self.setRowShow(self.tr('Additional user arrays'))
        self.setRowShow(self.tr('Volume solution control'), True)
        self.setRowShow(self.tr('Surface solution control'), True)
        self.setRowShow(self.tr('Lagrangian solution control'), m_lagr)
        self.setRowShow(self.tr('Profiles'), True)
        self.setRowShow(self.tr('Balance by zone'), (not is_ncfd))

        # Calculation management

        self.setRowShow(self.tr('Performance tuning'), True)
        self.setRowShow(self.tr('Prepare batch calculation'), True)

        # Update boundary zones display
        boundary_zone_labels = LocalizationModel("BoundaryZone", case).getSortedZoneLabels()
        self.updateBrowserZones(boundary_zone_labels, "Boundary conditions")

        # Update volume zones display
        volume_zone_labels = LocalizationModel("VolumicZone", case).getSortedZoneLabels()
        self.updateBrowserZones(volume_zone_labels, "Volume conditions")

        self.__hideRow()

    def __hideRow(self):
        """Only for developement purpose"""

    def tr(self, text):
        """
        Translation
        """
        return QCoreApplication.translate('BrowserView', text)


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
