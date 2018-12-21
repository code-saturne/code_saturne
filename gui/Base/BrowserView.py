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
from code_saturne.Base.QtWidgets import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.Base.BrowserForm import Ui_BrowserForm
from code_saturne.Base.Toolbox import GuiParam, displaySelectedPage
from code_saturne.Base.QtPage import to_qvariant, from_qvariant, to_text_string

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("BrowserView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

class TreeItem:
    def __init__(self, data, typename, parent=None):
        self.parentItem = parent
        self.itemData = data
        self.itemType = typename
        self.itemIcon = None
        self.childItems = []

    def appendChild(self, item):
        self.childItems.append(item)

    def child(self, row):
        return self.childItems[row]

    def childCount(self):
        return len(self.childItems)

    def columnCount(self):
        return len(self.itemData)

    def data(self, column):
        return self.itemData[column]

    def parent(self):
        return self.parentItem

    def row(self):
        if self.parentItem:
            return self.parentItem.childItems.index(self)

        return 0

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

        rootData = []
        rootData.append(to_qvariant("Pages"))

        self.rootItem = TreeItem(rootData, "folder")
        self.populateModel(data.split("\n"), self.rootItem)


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
            return to_qvariant()

        item = index.internalPointer()
        column = index.column()

        if role == Qt.DisplayRole:
            # return text for columns
            if column == 0:
                return to_qvariant(item.itemData[column])

        elif role == Qt.DecorationRole:
            # return icon for first column
            if column == 0:
                style = QWidget().style()
                if item.itemType == "folder-new":
                    icon = style.standardIcon(QStyle.SP_FileDialogNewFolder)
                elif item.itemType == "folder-close":
                    icon = style.standardIcon(QStyle.SP_DirClosedIcon)
                elif item.itemType == "folder-open":
                    icon = style.standardIcon(QStyle.SP_DirOpenIcon)
                elif item.itemType == "file-open":
                    icon = style.standardIcon(QStyle.SP_FileIcon)
                elif item.itemType == "file-new":
                    icon = style.standardIcon(QStyle.SP_FileLinkIcon)
                if sys.platform.startswith("win"):
                    if item.itemType == "file-open":
                        icon = style.standardIcon(QStyle.SP_ToolBarHorizontalExtensionButton)
                    elif item.itemType == "file-new":
                        icon = style.standardIcon(QStyle.SP_ToolBarHorizontalExtensionButton)
                return to_qvariant(icon)

        return to_qvariant()


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

        return to_qvariant()


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

        #FIXME: why childItem can be None?
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
                raise ValueError("This flags is not implemented")

            if self.hasChildren(index):
                result += self.match(self.index(0, index.column(), index),
                                     role, value, hits, flags)

        return result


    def itemLocalization(self, data, role=Qt.DisplayRole):
        """
        """
        info = []
        search_item = from_qvariant(to_qvariant(data), to_text_string)

        start = self.index(0, 0, QModelIndex())
        indexList = self.match(start, role, search_item, -1, Qt.MatchExactly)

        for index in indexList:
            item   = index.internalPointer()
            column = index.column()
            row    = index.row()
            parent = self.parent(index)

            info.append( (row, column, parent) )

        return info


    def populateModel(self, lines, parent):
        """
        @type lines: C{QString}
        @param lines:
        @type parent: C{QModelIndex}
        @param parent: parent of the item
        """
        parents = []
        indentations = []

        parents.append(parent)
        indentations.append(0)

        for number in range(len(lines)):
            position = 0
            while position < len(lines[number]):
                if lines[number][position] != " ":
                    break
                position += 1

            lineData = lines[number][position:].strip()

            if lineData:
                # Read the column data from the rest of the line.
                columnStrings = lineData.split("\t")

                columnData = []
                for column in range(0, len(columnStrings)):
                    columnData.append(columnStrings[column])

                if position == 0:
                    typename = "folder-new"
                else:
                    typename = "file-new"

                if position > indentations[-1]:
                    # The last child of the current parent is now the new parent
                    # unless the current parent has no children.
                    if parents[-1].childCount() > 0:
                        parents.append(parents[-1].child(parents[-1].childCount() - 1))
                        indentations.append(position)

                else:
                    while position < indentations[-1] and len(parents) > 0:
                        parents.pop()
                        indentations.pop()

                # Append a new item to the current parent's list of children.
                parents[-1].appendChild(TreeItem(columnData, typename, parents[-1]))

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

class BrowserView(QWidget, Ui_BrowserForm):
    """
    Class for the browser widget
    """
    def __init__(self):
        """
        Constructor
        """
        QWidget.__init__(self)

        Ui_BrowserForm.__init__(self)
        self.setupUi(self)

        tree = self._browser()
        self.model = TreeModel(from_qvariant(to_qvariant(tree), to_text_string))

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


    def _browser(self):
        tree ="""
    Identity and paths
Calculation environment
    Meshes selection
    Notebook
Thermophysical models
    Calculation features
    Main fields
    Deformable mesh
    Turbulence models
    Thermal model
    Gas combustion
    Pulverized fuel combustion
    Electrical models
    Radiative transfers
    Conjugate heat transfer
    Atmospheric flows
    Species transport
    Turbomachinery
    Groundwater flows
    Fans
    Non condensable gases
    Thermodynamics
Closure modeling
    Interfacial momentum transfer
    Interfacial area
    Interfacial enthalpy transfer
    Nucleate boiling parameters
    Droplet condensation-evaporation
    Particles interactions
Physical properties
    Reference values
    Fluid properties
    Gravity
Particles and droplets tracking
    Global settings
    Statistics
Volume conditions
    Volume regions definition
    Main fields initialization
    Initialization
    Head losses
    Porosity
    Source terms
    Coriolis Source Terms
    Groundwater laws
Boundary conditions
    Boundary regions definition
    Boundary conditions
    Particle boundary conditions
    Fluid structure interaction
    Cathare Coupling
Numerical parameters
    Global parameters
    Equation parameters
    Time settings
Calculation control
    Time averages
    Additional user arrays
    Output control
    Volume solution control
    Surface solution control
    Lagrangian solution control
    Profiles
    Balance by zone
Calculation management
    Start/Restart
    Performance tuning
    OpenTurns study
    Prepare batch calculation
"""
        return tree


    def setRowClose(self, string):
        log.debug("setRowClose(): %s" % string)
        itemInfoList = self.model.itemLocalization(string)
        for itemInfo in itemInfoList:
            row    = itemInfo[0]
            column = itemInfo[1]
            parent = itemInfo[2]
            self.treeView.setRowHidden(row, parent, True)


    def setRowOpen(self, string):
        log.debug("setRowOpen(): %s" % string)
        itemInfoList = self.model.itemLocalization(string)
        for itemInfo in itemInfoList:
            row    = itemInfo[0]
            column = itemInfo[1]
            parent = itemInfo[2]
            self.treeView.setRowHidden(row, parent, False)


    def setRowShow(self, string, status=True):
        log.debug("setRowVisible(): %s" % string)
        itemInfoList = self.model.itemLocalization(string)
        hidden = not status
        for itemInfo in itemInfoList:
            row    = itemInfo[0]
            column = itemInfo[1]
            parent = itemInfo[2]
            self.treeView.setRowHidden(row, parent, hidden)


    def isRowClose(self, string):
        log.debug("isRowClose(): %s" % string)
        itemInfoList = self.model.itemLocalization(string)
        for itemInfo in itemInfoList:
            row    = itemInfo[0]
            column = itemInfo[1]
            parent = itemInfo[2]
            index  = self.model.index(row, column, parent)
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
        #self.actionExpand.setShortcut(self.tr("F5"))
        self.actionExpand.triggered.connect(self.openTreeFolder)

        self.actionCollapse = QAction(self.tr("Collapse"), self.treeView)
        #self.actionCollapse.setShortcut(self.tr("F6"))
        self.actionCollapse.triggered.connect(self.closeTreeFolder)

        # ... TODO
        #self.actionWelcome = QAction(self.tr("Welcome page"), self.treeView)

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


    def display(self, root, case, stbar, study, tree):
        """
        """
        index = self.treeView.currentIndex()
        item  = index.internalPointer()
        name  = item.itemData[0]
        case['current_tab'] = 0
        case['current_index'] = index
        return displaySelectedPage(name, root, case, stbar, study, tree)


    def isFolder(self):
        """
        Return True if current item is a folder (parent)
        """
        index = self.treeView.currentIndex()
        item  = index.internalPointer()
        return item.childCount() != 0


    def openSingleFolder(self, string):
        """
        Open a single folder of the Tree.
        """
        itemInfoList = self.model.itemLocalization(string)
        for itemInfo in itemInfoList:
            row    = itemInfo[0]
            column = itemInfo[1]
            parent = itemInfo[2]
            index  = self.model.index(row, column, parent)
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
            row    = itemInfo[0]
            column = itemInfo[1]
            parent = itemInfo[2]
            index  = self.model.index(row, column, parent)
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

        self.setRowClose(self.tr('Notebook'))

        self.setRowClose(self.tr('Thermophysical models'))
        """
        self.setRowClose(self.tr('Calculation features'))
        self.setRowClose(self.tr('Main fields'))
        self.setRowClose(self.tr('Deformable mesh'))
        self.setRowClose(self.tr('Turbulence models'))
        self.setRowClose(self.tr('Thermal model'))
        self.setRowClose(self.tr('Gas combustion'))
        self.setRowClose(self.tr('Pulverized fuel combustion'))
        self.setRowClose(self.tr('Electrical models'))
        self.setRowClose(self.tr('Radiative transfers'))
        self.setRowClose(self.tr('Conjugate heat transfer'))
        self.setRowClose(self.tr('Atmospheric flows'))
        self.setRowClose(self.tr('Species transport'))
        self.setRowClose(self.tr('Turbomachinery'))
        self.setRowClose(self.tr('Groundwater flows'))
        self.setRowClose(self.tr('Fans'))
        self.setRowClose(self.tr('Non condensable gases'))
        self.setRowClose(self.tr('Thermodynamics'))
        """

        self.setRowClose(self.tr('Closure modeling'))
        """
        self.setRowClose(self.tr('Interfacial momentum transfer'))
        self.setRowClose(self.tr('Interfacial area'))
        self.setRowClose(self.tr('Interfacial enthalpy transfer'))
        self.setRowClose(self.tr('Nucleate boiling parameters'))
        self.setRowClose(self.tr('Droplet condensation-evaporation'))
        self.setRowClose(self.tr('Particles interactions'))
        """

        self.setRowClose(self.tr('Physical properties'))
        """
        self.setRowClose(self.tr('Reference values'))
        self.setRowClose(self.tr('Fluid properties'))
        self.setRowClose(self.tr('Gravity'))
        """

        self.setRowClose(self.tr('Particles and droplets tracking'))
        """
        self.setRowClose(self.tr('Global settings'))
        self.setRowClose(self.tr('Statistics'))
        """

        self.setRowShow(self.tr('Volume conditions'), False)
        """
        self.setRowClose(self.tr('Volume regions definition'))
        self.setRowClose(self.tr('Main fields initialization'))
        self.setRowClose(self.tr('Initialization'))
        self.setRowClose(self.tr('Head losses'))
        self.setRowClose(self.tr('Porosity'))
        self.setRowClose(self.tr('Source terms'))
        self.setRowClose(self.tr('Coriolis Source Terms'))
        self.setRowClose(self.tr('Groundwater laws'))
        """

        self.setRowClose(self.tr('Boundary conditions'))
        self.setRowClose(self.tr('Particle boundary conditions'))
        self.setRowClose(self.tr('Fluid structure interaction'))
        self.setRowClose(self.tr('Cathare Coupling'))

        self.setRowClose(self.tr('Numerical parameters'))
        """
        self.setRowClose(self.tr('Global parameters'))
        self.setRowClose(self.tr('Equation parameters'))
        self.setRowClose(self.tr('Time settings'))
        """

        self.setRowClose(self.tr('Time averages'))
        self.setRowClose(self.tr('Additional user arrays'))
        self.setRowShow(self.tr('Output control'))
        self.setRowClose(self.tr('Volume solution control'))
        self.setRowClose(self.tr('Surface solution control'))
        self.setRowClose(self.tr('Lagrangian solution control'))
        self.setRowClose(self.tr('Profiles'))
        self.setRowClose(self.tr('Balance by zone'))

        self.setRowClose(self.tr('Start/Restart'))
        self.setRowClose(self.tr('OpenTurns study'))

        self.__hideRow()


    def configureTree(self, case):
        """
        Public method.
        Configures the browser with users data.
        """

        if case['prepro'] == True:
            return self.__configureTreePrepro(case)

        p_module = ''
        if case.xmlRootNode().tagName == "NEPTUNE_CFD_GUI":
            p_module = 'neptune_cfd'

        # Precompute some values
        #-----------------------

        m_tbm = False
        m_fans = False
        m_lagr = False
        m_ale = 0
        m_thermal = 0
        m_cht = False
        m_gas_comb = False
        m_sf_comb = False
        m_elec = False
        m_atmo = False
        m_rad = False
        m_comp = False
        m_gwf = False

        node_pm = case.xmlGetNode('thermophysical_models')

        if node_pm:

            node = node_pm.xmlGetNode('ale_method', 'status')
            if node and node['status'] == 'on':
                m_ale = 1
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
                node = node_pm.xmlGetNode('atmospheric_flows',  'model')
                if node and node['model'] != 'off':
                    m_atmo = True
                    if node['model'] == 'constant':
                        m_thermal = -1
                    elif node['model'] in ('dry', 'humid'):
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

            node = node_pm.xmlGetNode('groundwater_model',  'model')
            if node and node['model'] != 'off':
                m_gwf = True
                m_thermal = -1
                m_ale = -1
                m_fans = False

            if m_thermal > 0:
                m_rad = True
                m_cht = True

        node = case.xmlGetNode('lagrangian', 'model')
        if node and node['model'] != "off":
            m_lagr = True

        # Options for NEPTUNE_CFD if present:

        m_ncfd = {}
        m_ncfd['non_condens'] = False
        m_ncfd['nucleate_boiling'] = False
        m_ncfd['droplet_condens'] = False
        m_ncfd['particles_interactions'] = False
        m_ncfd['itf_area'] = False
        m_ncfd['itf_h_transfer'] = False

        ncfd_fields = 0

        if p_module == 'neptune_cfd' and node_pm:

            m_ale = -1
            m_fans = False
            m_thermal = -1
            m_rad = True
            m_cht = True

            node_f = node_pm.xmlInitNode('fields')
            fields = node_f.xmlGetNodeList('field')
            ncfd_fields = len(fields)

        if ncfd_fields > 1:
            from code_saturne.Pages.MainFieldsModel import MainFieldsModel
            from code_saturne.Pages.NonCondensableModel import NonCondensableModel
            from code_saturne.Pages.InterfacialForcesModel import InterfacialForcesModel
            predefined_flow = MainFieldsModel(case).getPredefinedFlow()

            if (len(MainFieldsModel(case).getSolidFieldIdList()) > 0):
                m_ncfd['particles_interactions'] = True
            if (len(MainFieldsModel(case).getDispersedFieldList()) > 0
                or InterfacialForcesModel(case).getBubblesForLIMStatus() == 'on'):
                m_ncfd['itf_area'] = True

            if predefined_flow == "free_surface":
                m_ncfd['non_condens'] = True
                m_ncfd['nucleate_boiling'] = True
                m_ncfd['itf_h_transfer'] = True
            elif predefined_flow == "boiling_flow":
                m_ncfd['non_condens'] = True
                m_ncfd['nucleate_boiling'] = True
                m_ncfd['itf_h_transfer'] = True
            elif predefined_flow == "droplet_flow":
                m_ncfd['non_condens'] = True
                m_ncfd['droplet_condens'] = True
                m_ncfd['itf_h_transfer'] = True
            elif predefined_flow == "particles_flow":
                m_ncfd['particles_interactions'] = True
                m_ncfd['itf_h_transfer'] = True

        is_ncfd = (p_module == 'neptune_cfd')

        # Manage visibility
        #------------------

        self.setRowOpen(self.tr('Notebook'))

        # Thermophysical Models

        self.setRowShow(self.tr('Thermophysical models'), True)
        self.setRowShow(self.tr('Calculation features'), True)
        self.setRowShow(self.tr('Main fields'), (p_module == 'neptune_cfd'))
        self.setRowShow(self.tr('Deformable mesh'), (m_ale > -1))
        self.setRowShow(self.tr('Turbulence models'))
        self.setRowShow(self.tr('Thermal model'), (m_thermal > -1))
        self.setRowShow(self.tr('Gas combustion'), m_gas_comb)
        self.setRowShow(self.tr('Pulverized fuel combustion'), m_sf_comb)
        self.setRowShow(self.tr('Electrical models'), m_elec)
        self.setRowShow(self.tr('Radiative transfers'), m_rad)
        self.setRowShow(self.tr('Conjugate heat transfer'), m_cht)
        self.setRowShow(self.tr('Atmospheric flows'), m_atmo)
        self.setRowShow(self.tr('Species transport'))
        self.setRowShow(self.tr('Turbomachinery'), m_tbm)
        self.setRowShow(self.tr('Groundwater flows'), m_gwf)
        self.setRowShow(self.tr('Fans'), m_fans)

        self.setRowShow(self.tr('Non condensable gases'), m_ncfd['non_condens'])
        self.setRowShow(self.tr('Thermodynamics'), is_ncfd)

        # Closure modeling

        self.setRowShow(self.tr('Closure modeling'), (ncfd_fields > 1))
        self.setRowShow(self.tr('Interfacial momentum transfer'), (ncfd_fields > 1))
        self.setRowShow(self.tr('Interfacial enthalpy transfer'), m_ncfd['itf_h_transfer'])
        self.setRowShow(self.tr('Interfacial area'), m_ncfd['itf_area'])
        self.setRowShow(self.tr('Nucleate boiling parameters'), m_ncfd['nucleate_boiling'])
        self.setRowShow(self.tr('Droplet condensation-evaporation'), m_ncfd['droplet_condens'])
        self.setRowShow(self.tr('Particles interactions'), m_ncfd['particles_interactions'])

        # Physical properties

        self.setRowShow(self.tr('Physical properties'), (not m_gwf))
        self.setRowShow(self.tr('Reference values'), (not (m_gwf or is_ncfd)))
        self.setRowShow(self.tr('Fluid properties'), (not (m_gwf or is_ncfd)))
        self.setRowShow(self.tr('Gravity'), (not m_gwf))

        # Particles and droplets tracking

        self.setRowShow(self.tr('Particles and droplets tracking'), m_lagr)
        self.setRowShow(self.tr('Global settings'), m_lagr)
        self.setRowShow(self.tr('Statistics'), m_lagr)

        # Volume conditions

        self.setRowShow(self.tr('Volume conditions'), True)
        self.setRowShow(self.tr('Volume regions definition'), True)

        node_domain = case.xmlGetNode('solution_domain')
        node_vol = node_domain.xmlGetNode('volumic_conditions')
        init = False
        z_st = False
        z_head_loss = False
        z_porosity = False
        z_groundwater = False

        for node in node_vol.xmlGetChildNodeList('zone'):
            if (node['initialization'] == 'on'):
                init = True
            if (node['momentum_source_term'] == 'on'
                or node['mass_source_term'] == 'on'
                or node['thermal_source_term'] == 'on'
                or node['scalar_source_term'] == 'on'):
                z_st = True
            if node['head_losses'] == 'on':
                z_head_loss = True
            if node['porosity'] == 'on':
                z_porosity = True
            if node['groundwater_law'] == 'on':
                z_groundwater = True

        self.setRowShow(self.tr('Main fields initialization'), is_ncfd and init)
        self.setRowShow(self.tr('Initialization'), (not is_ncfd) and init)
        self.setRowShow(self.tr('Head losses'), z_head_loss)
        self.setRowShow(self.tr('Porosity'), z_porosity)
        self.setRowShow(self.tr('Source terms'), z_st)
        self.setRowShow(self.tr('Coriolis Source Terms'), (not (m_gwf or is_ncfd)))
        self.setRowShow(self.tr('Groundwater laws'), z_groundwater)

        # Boundary conditions

        self.setRowShow(self.tr('Boundary conditions'))
        self.setRowShow(self.tr('Boundary regions definition'), True)
        self.setRowShow(self.tr('Particle boundary conditions'), m_lagr)
        self.setRowShow(self.tr('Fluid structure interaction'), (m_ale > 0))
        self.setRowShow(self.tr('Cathare Coupling'), is_ncfd)

        # Numerical parameters

        self.setRowShow(self.tr('Numerical parameters'))
        self.setRowShow(self.tr('Global parameters'))
        self.setRowShow(self.tr('Equation parameters'))
        self.setRowShow(self.tr('Time settings'))

        # Calculation control

        self.setRowShow(self.tr('Time averages'), True)
        self.setRowShow(self.tr('Additional user arrays'))
        self.setRowShow(self.tr('Output control'), True)
        self.setRowShow(self.tr('Volume solution control'), True)
        self.setRowShow(self.tr('Surface solution control'), True)
        self.setRowShow(self.tr('Lagrangian solution control'), m_lagr)
        self.setRowShow(self.tr('Profiles'), True)
        self.setRowShow(self.tr('Balance by zone'), (not is_ncfd))

        # Calculation management

        self.setRowShow(self.tr('Start/Restart'), True)
        self.setRowShow(self.tr('Performance tuning'), True)
        self.setRowShow(self.tr('Prepare batch calculation'), True)
        self.setRowShow(self.tr('OpenTurns study'), case['salome'])

        # End of test of physical module

        self.__hideRow()


    def __hideRow(self):
        """Only for developement purpose"""


    def tr(self, text):
        """
        Translation
        """
        return text

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
