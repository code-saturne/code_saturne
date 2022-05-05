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
This module defines basic classes used for the Pages construction.

This module defines the following classes:
- ComboModel
- IntValidator
- DoubleValidator
- RegExpValidator
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import sys
import os
import logging
import locale

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

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("QtPage")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Compatibility PyQt API #1 and #2 and PySide
#-------------------------------------------------------------------------------

#==============================================================================
# Data types
#==============================================================================

TEXT_TYPES = (str,)

#==============================================================================
# Classes used to handle ComboBox subsections (groups)
#==============================================================================

# Qt Role used to differentiate whether a ComboBox entry is a group child
# or not
GroupRole = Qt.UserRole


class GroupDelegate(QStyledItemDelegate):
    """
    Class used to create a group delegate which applies a tab to a combobox
    entry if it is a group child. Only applied to the visual style!
    """
    def initStyleOption(self, option, index):
        super(GroupDelegate, self).initStyleOption(option, index)
        if not index.data(GroupRole):
            option.text = "   " + option.text

class GroupItem():
    """
    GroupItem class.
    Handles the subsections in a ComboModel.
    """

    def __init__(self, name, index):
        """
        Init method. Inputs are:
        name  : String with the name of the group
        index : Position of the group name in the combobox list of elements
        """
        self._name  = name
        self._index = index

        self._children = []
        self._number_of_children = 0


    def addChild(self, name):
        """
        Method to add a child with the given name.
        Returns the index of the child to be used if the 'insertRow' method.
        """
        self._children.append(name)
        self._number_of_children += 1

        child_index = self._index + self._number_of_children

        return child_index

    def getIndex(self):
        """
        Method which returns the index of the group
        """
        return self._index

    def getChildIndex(self, childName):
        """
        Method to retrieve the index of a given child
        """
        if childName not in self._children:
            msg = "%s is not part of group %s\n" % (childName, self._name)
            raise Exception(msg)

        child_index = self._children.index(childName) + self._index + 1
        return child_index


    def getNumberOfChildren(self):
        """
        Method which returns the number of children in the group
        """
        return self._number_of_children


    def incrementIndex(self):

        self._index += 1

    def generateItem(self):

        item  = QStandardItem(str(self._name))
        item.setData(True, GroupRole)

        font = item.font()
        font.setBold(True)
        font.setItalic(True)
        item.setFont(font)
        item.setFlags(item.flags() & ~Qt.ItemIsSelectable)

        return item
#==============================================================================
# Strings
#==============================================================================

def is_text_string(obj):
    """Return True if `obj` is a text string,
              False if it is anything else,
                    like binary data (Python 3) or
                    QString (Python 2, PyQt API #1)"""
    return isinstance(obj, str)

def to_text_string(obj, encoding=None):
    """Convert `obj` to (unicode) text string"""
    if encoding is None:
        return str(obj)
    elif isinstance(obj, str):
        # In case this function is not used properly, this could happen
        return obj
    else:
        return str(obj, encoding)

#==============================================================================
# QVariant conversion utilities
#
# Note: to_qvariant was removed recently; from_qvariant may be removed in
#       many_places and the "to_text_string" variants could be replaced
#       by "to_str"; where the value is already known to be of the correct
#       type, or None handled in the associated code,
#       this could even be ignored
#==============================================================================

import collections

if os.environ.get('QT_API', 'pyqt') == 'pyqt':

    def from_qvariant(qobj=None, convfunc=None):
        """Convert QVariant object to Python object
        This is a transitional function from PyQt API #1 (QVariant exists)
        to PyQt API #2 and Pyside (QVariant does not exist)"""
        if (qobj != None):
            if convfunc in TEXT_TYPES or convfunc is to_text_string:
                return str(qobj)
            elif convfunc is int:
                return int(qobj)
            elif convfunc is float:
                try:
                    return float(qobj)
                except Exception:
                    return locale.atof(qobj)
            else:
                return qobj
        else:
            return qobj

else:

    def from_qvariant(qobj=None, pytype=None):
        """Convert QVariant object to Python object
        This is a transitional function from PyQt API #1 (QVariant exist)
        to PyQt API #2 and Pyside (QVariant does not exist)"""
        return qobj

def qbytearray_to_str(qba):
    """Convert QByteArray object to str in a way compatible with Python 2/3"""
    return str(bytes(qba.toHex().data()).decode())


def to_str(s):
    """Convert argument to string, using an empty string for None"""
    if s is None:
        return ''
    return str(s)

#==============================================================================
# Wrappers around QFileDialog static methods
#==============================================================================

def getexistingdirectory(parent=None, caption='', basedir='',
                         options=QFileDialog.ShowDirsOnly):
    """Wrapper around QtGui.QFileDialog.getExistingDirectory static method
    Compatible with PyQt >=v4.4 (API #1 and #2) and PySide >=v1.0"""
    # Calling QFileDialog static method
    if sys.platform == "win32":
        # On Windows platforms: redirect standard outputs
        _temp1, _temp2 = sys.stdout, sys.stderr
        sys.stdout, sys.stderr = None, None
    try:
        result = QFileDialog.getExistingDirectory(parent, caption, basedir,
                                                  options)
    finally:
        if sys.platform == "win32":
            # On Windows platforms: restore standard outputs
            sys.stdout, sys.stderr = _temp1, _temp2
    if not is_text_string(result):
        # PyQt API #1
        result = to_text_string(result)
    return result


def _qfiledialog_wrapper(attr, parent=None, caption='', basedir='',
                         filters='', selectedfilter='', options=None):
    if options is None:
        options = QFileDialog.Options(0)

    try:
        from code_saturne.gui.base.QtCore import QString
    except ImportError:
        QString = None  # analysis:ignore

    tuple_returned = True
    try:
        func = getattr(QFileDialog, attr+'AndFilter')
    except AttributeError:
        func = getattr(QFileDialog, attr)
        if QString is not None:
            selectedfilter = QString()
            tuple_returned = False

    if sys.platform == "win32":
        # On Windows platforms: redirect standard outputs
        _temp1, _temp2 = sys.stdout, sys.stderr
        sys.stdout, sys.stderr = None, None
    try:
        result = func(parent, caption, basedir,
                      filters, selectedfilter, options)
    except TypeError:
        result = func(parent, caption, basedir, filters, options)
    finally:
        if sys.platform == "win32":
            # On Windows platforms: restore standard outputs
            sys.stdout, sys.stderr = _temp1, _temp2

    # Processing output
    if tuple_returned:
        output, selectedfilter = result
    else:
        output = result
    if QString is not None:
        # PyQt API #1: conversions needed from QString/QStringList
        selectedfilter = to_text_string(selectedfilter)
        if isinstance(output, QString):
            # Single filename
            output = to_text_string(output)
        else:
            # List of filenames
            output = [to_text_string(fname) for fname in output]

    # Always returns the tuple (output, selectedfilter)
    return output, selectedfilter


def getopenfilename(parent=None,
                    caption='',
                    basedir='',
                    filters='',
                    selectedfilter='',
                    options=None):
    return _qfiledialog_wrapper('getOpenFileName',
                                parent=parent,
                                caption=caption,
                                basedir=basedir,
                                filters=filters,
                                selectedfilter=selectedfilter,
                                options=options)


def getopenfilenames(parent=None,
                     caption='',
                     basedir='',
                     filters='',
                     selectedfilter='',
                     options=None):
    return _qfiledialog_wrapper('getOpenFileNames',
                                parent=parent,
                                caption=caption,
                                basedir=basedir,
                                filters=filters,
                                selectedfilter=selectedfilter,
                                options=options)


def getsavefilename(parent=None,
                    caption='',
                    basedir='',
                    filters='',
                    selectedfilter='',
                    options=None):
    return _qfiledialog_wrapper('getSaveFileName',
                                parent=parent,
                                caption=caption,
                                basedir=basedir,
                                filters=filters,
                                selectedfilter=selectedfilter,
                                options=options)


#-------------------------------------------------------------------------------
# QComboBox model
#-------------------------------------------------------------------------------

class ComboModel:
    """
    Class to build a model (QStandardItemModel) used with a QComboBox.

    Main attributes of class are:

    combo: QComboBox passed as arguments. It uses the model
    model: QStandardItemModel which contains items

    dicoV2M: correspondance between strings in Qt view and strings in parameters
    dicoM2V: correspondance between strings in parameters and strings in Qt view

    items: tuple which contains all model strings (usefull to get its index in the model)
    """
    def __init__(self, combo, rows=0, columns=0):
        """
        Initialization
        """
        self.combo   = combo

        self.rows    = rows
        self.columns = columns
        self.last    = 0

        self.model   = QStandardItemModel()
        self.model.clear()
        self.model.setRowCount(rows)
        self.model.setColumnCount(columns)

        self.dicoV2M = {}
        self.dicoM2V = {}

        self.item_groups = {}

        self.items   = []
        self.combo.setModel(self.model)
        self.combo.setItemDelegate(GroupDelegate(self.combo))


    def addItemGroup(self, group_name):

        if group_name in self.item_groups.keys():
            return

        index = self.last
        self.item_groups[group_name] = GroupItem(group_name, index)

        item = self.item_groups[group_name].generateItem()
        self.model.setItem(index, item)

        self.items.append(group_name)

        self.last += 1

    def addItemList(self, list_of_views_and_models, warn=False, groupName=None):
        """
        Insert of a list of items in the model.
        list_of_views_and_models is a list of 2-element list
          e.g. [ ["View1", "Model1"], ["View2", "Model2"] ]
        """
        for view, model in list_of_views_and_models:
            self.addItem(view, model, warn, groupName)

    def hasItem(self, index=None, str_view="", str_model=""):
        if index is not None:
            return index < self.last

        elif str_model:
            return (str_model in self.items)

        elif str_view:
            return (str_view in self.dicoV2M.keys())

    def addItem(self, str_view, str_model="", warn=False, groupName=None):
        """
        Insert an item in the model.

        str_view: string to be displayed in the view.
        For example, 'Eulerian/Lagrangian Multi-phase Treatment'

        str_model: correponding string used in the model.
        For example, 'lagrangian'

        warn: If True, entry is marked with a color.
        """
        if self.hasItem(str_view=str_view):
            return

        item = QStandardItem(str(str_view))
        item.setData(True, GroupRole)

        if groupName in self.item_groups.keys():
            # Insert the child after the group name in the ComboModel
            item.setData(False, GroupRole)
            index = self.item_groups[groupName].addChild(str(str_view))
            self.model.setItem(index, item)

            # In case groups were created before, update the index of
            # the following groups in the list
            for key in self.item_groups.keys():
                gpe = self.item_groups[key]
                if gpe.getIndex() >= index:
                    gpe.incrementIndex()

                    idx = gpe.getIndex()
                    it  = gpe.generateItem()
                    self.model.setItem(idx, it)

        else:
            index = self.last
            self.model.setItem(index, item)

        if warn:
            self.combo.setItemData(index,
                                   QColor(Qt.red),
                                   Qt.TextColorRole)

        self.last += 1

        if not str_model: str_model = str_view

        self.items.append(str_model)

        self.dicoM2V[str_model] = str_view
        self.dicoV2M[str_view]  = str_model


    def modifyItem(self, old_str_view, new_str_view, new_str_model=""):
        """
        Modify string names.
        """
        if old_str_view in self.items:
            index = self.items.index(old_str_view)
            self.items[index] = new_str_view

            old_str_model = dicoV2M[old_str_view]
            if new_str_model == "":
                new_str_model = old_str_model

            del self.dicoM2V[str_model]
            del self.dicoV2M[str_view]

            self.dicoM2V[new_str_model] = new_str_view
            self.dicoV2M[new_str_view] = new_str_model

    def flushItems(self):
        """
        Remove all items
        """
        while self.last > 0:
            self.delItem(0)

    def delItem(self, index=None, str_model="", str_view=""):
        """
        Remove the item specified with its index or a string.
        """
        if index is not None:
            self.__deleteItem(index)

        elif str_model:
            index = self.items.index(str_model)
            self.__deleteItem(index)

        elif str_view:
            str_model = self.dicoV2M[str_view]
            index = self.items.index(str_model)
            self.__deleteItem(index)
        self.combo.removeItem(index)
        self.last = self.last - 1


    def __deleteItem(self, index):
        """
        Delete the item specified with its index
        """
        str_model  = self.items[index]
        str_view = self.dicoM2V[str_model]
        del self.items[index]
        del self.dicoV2M[str_view]
        del self.dicoM2V[str_model]


    def __disableItem(self, index):
        """
        Disable the item specified with its index
        """
        self.model.item(index).setEnabled(False)


    def __enableItem(self, index):
        """
        Enable the item specified with its index
        """
        self.model.item(index).setEnabled(True)


    def disableItem(self, index=None, str_model="", str_view=""):
        """
        Disable the item specified with its index or a string.
        """
        if index is not None:
            self.__disableItem(index)

        elif str_model:
            index = self.items.index(str_model)
            self.__disableItem(index)

        elif str_view:
            str_model = self.dicoV2M[str_view]
            index = self.items.index(str_model)
            self.__disableItem(index)


    def enableItem(self, index=None, str_model="", str_view=""):
        """
        Enable the item specified with its index or a string.
        """
        if index is not None:
            self.__enableItem(index)

        elif str_model:
            index = self.items.index(str_model)
            self.__enableItem(index)

        elif str_view:
            str_model = self.dicoV2M[str_view]
            index = self.items.index(str_model)
            self.__enableItem(index)


    def setItem(self, index=None, str_model="", str_view=""):
        """
        Set item as current.
        """
        if index is not None:
            self.combo.setCurrentIndex(index)

        elif str_model:
            try:
                index = self.items.index(str_model)
                self.combo.setCurrentIndex(index)
            except Exception:
                self._displayWarning(str_model)
                # Throw signals to ensure XML is updated (not very elegant)
                self.combo.activated[str].emit(self.dicoM2V[self.items[0]])
                self.combo.currentTextChanged[str].emit(self.dicoM2V[self.items[0]])
                self.combo.currentIndexChanged[int].emit(0)
                index = 0

        elif str_view:
            try:
                str_model = self.dicoV2M[str_view]
                index = self.items.index(str_model)
                self.combo.setCurrentIndex(index)
            except Exception:
                self._displayWarning(str_model)
                index = 0

    def _displayWarning(self, str_model):
        if self.combo.accessibleName() != "":
            name = self.combo.accessibleName()
        else:
            name = self.combo.objectName()
        title = "Warning in " + name
        msg = str_model + " is not in list: " + str(self.items) + "\n"
        msg += "Value reset to " + str(self.items[0]) + " but should be checked.\n"
        QMessageBox.warning(self.combo, title, msg)

    def enableItem(self, index=None, str_model="", str_view=""):
        """
        Enable the item specified with its index or a string.
        """
        if index is not None:
            self.__enableItem(index)

        elif str_model:
            index = self.items.index(str_model)
            self.__enableItem(index)

        elif str_view:
            str_model = self.dicoV2M[str_view]
            index = self.items.index(str_model)
            self.__enableItem(index)


    def hide(self):
        """
        Hide combobox.
        """
        self.combo.hide()


    def show(self):
        """
        Show combobox.
        """
        self.combo.show()


    def getItems(self):
        """
        Get the tuple of items.
        """
        return self.items

#-------------------------------------------------------------------------------
# Validators for editors
#-------------------------------------------------------------------------------

vmax = 2147483647
vmin = -vmax

class IntValidator(QIntValidator):
    """
    Validator for integer data.
    """
    def __init__(self, parent, min=vmin, max=vmax):
        """
        Initialization for validator
        """
        QIntValidator.__init__(self, parent)
        self.parent = parent
        self.state = QValidator.Invalid
        self.__min = min
        self.__max = max

        if type(min) != int or type(max) != int:
            raise ValueError("The given parameters are not integers (warning: long are not allowed).")
        self.setBottom(min)
        self.setTop(max)

        self.exclusiveMin = False
        self.exclusiveMax = False
        self.exclusiveValues = []

        self.default = 0
        self.fix = False

        msg = ""
        if min > vmin and max == vmax:
            msg = self.tr("The integer value must be greater than or equal to %i" % min)
        elif min == vmin and max < vmax:
            msg = self.tr("The integer value must be lower than or equal to %i" % max)
        elif min > vmin and max < vmax:
            msg = self.tr("The integer value must be between %i and %i" % (min, max))

        self.parent.setStatusTip(msg)


    def setExclusiveMin(self, b=True):
        if type(b) != bool:
            raise ValueError("The given parameter is not a boolean.")
        self.exclusiveMin = b

        msg = ""
        if self.__min > vmin and self.__max == vmax:
            msg = self.tr("The integer value must be greater than %i" % self.__min)
        elif self.__min == vmin and self.__max < vmax:
            msg = self.tr("The integer value must be lower than or equal to %i" % self.__max)
        elif self.__min > vmin and self.__max < vmax:
            msg = self.tr("The integer value must be greater %i and lower than or equal to %i" % (self.__min, self.__max))

        self.parent.setStatusTip(msg)


    def setExclusiveMax(self, b=True):
        if type(b) != bool:
            raise ValueError("The given parameter is not a boolean.")
        self.exclusiveMax = b

        msg = ""
        if self.__min > vmin and self.__max == vmax:
            msg = self.tr("The integer value must be greater than or equal to %i" % self.__min)
        elif self.__min == vmin and self.__max < vmax:
            msg = self.tr("The integer value must be lower than %i" % self.__max)
        elif self.__min > vmin and self.__max < vmax:
            msg = self.tr("The integer value must be greater than or equal to %i and lower than %i" % (self.__min, self.__max))

        self.parent.setStatusTip(msg)


    def setExclusiveValues(self, l):
        if type(l) != list and type(l) != tuple:
            raise ValueError("The given parameter is not a list or a tuple.")
        self.exclusiveValues = l

        msg = ""
        for v in l:
            if self.__min > vmin or self.__max < vmax:
                msg = self.tr("All integers value must be greater than %i and lower than %i" % (self.__min, self.__max))

        self.parent.setStatusTip(msg)


    def setFixup(self, v):
        if type(v) != int:
            raise ValueError("The given parameter is not an integer.")
        self.default = v
        self.fix = True


    def fixup(self, stri):
        if self.fix:
            if not stri:
                stri = str(self.default)


    def validate(self, stri, pos):
        """
        Validation method.

        QValidator.Invalid       0  The string is clearly invalid.
        QValidator.Intermediate  1  The string is a plausible intermediate value during editing.
        QValidator.Acceptable    2  The string is acceptable as a final result; i.e. it is valid.
        """
        state = QIntValidator.validate(self, stri, pos)[0]

        try:
            x = from_qvariant(stri, int)
            valid = True
            pass
        except (TypeError, ValueError):
            x = 0
            valid = False

        if state == QValidator.Acceptable:
            if self.exclusiveMin and x == self.bottom():
                state = QValidator.Intermediate
            elif self.exclusiveMax and x == self.top():
                state = QValidator.Intermediate
            elif x in self.exclusiveValues:
                state = QValidator.Intermediate

        palette = self.parent.palette()

        if not valid or state == QValidator.Intermediate:
            palette.setColor(QPalette.Text, QColor("red"))
            self.parent.setPalette(palette)
        else:
            palette.setColor(QPalette.Text, QColor("black"))
            self.parent.setPalette(palette)

        self.state = state

        return (state, stri, pos)


    def tr(self, text):
        """
        """
        return text


class DoubleValidator(QDoubleValidator):
    """
    Validator for real data.
    """
    def __init__(self, parent, min=-1.e99, max=1.e99):
        """
        Initialization for validator
        """
        QDoubleValidator.__init__(self, parent)
        self.setLocale(QLocale(QLocale.C, QLocale.AnyCountry))
        self.parent = parent
        self.state = QValidator.Invalid
        self.__min = min
        self.__max = max

        self.setNotation(self.ScientificNotation)

        if type(min) != float or type(max) != float:
            raise ValueError("The given parameters are not floats.")
        self.setBottom(min)
        self.setTop(max)

        self.exclusiveMin = False
        self.exclusiveMax = False

        self.default = 0.0
        self.fix = False

        msg = ""
        if min > -1.e99 and max == 1.e99:
            msg = self.tr("The float value must be greater than %.2f" % min)
        elif min == -1.e99 and max < 1.e99:
            msg = self.tr("The float value must be lower than %.2f" % max)
        elif min > -1.e99 and max < 1.e99:
            msg = self.tr("The float value must be between %.2f and %.2f" % (min, max))

        self.parent.setStatusTip(str(msg))


    def setExclusiveMin(self, b=True):
        if type(b) != bool:
            raise ValueError("The given parameter is not a boolean.")
        self.exclusiveMin = b

        msg = ""
        if self.__min > -1.e99 and self.__max == 1.e99:
            msg = self.tr("The float value must be greater than %.1f" % self.__min)
        elif self.__min == -1.e99 and self.__max < 1.e99:
            msg = self.tr("The float value must be lower than or equal to %.1f" % self.__max)
        elif self.__min > -1.e99 and self.__max < 1.e99:
            msg = self.tr("The float value must be greater than %.1f and lower than or equal to %.1f" % (self.__min, self.__max))

        self.parent.setStatusTip(str(msg))


    def setExclusiveMax(self, b=True):
        if type(b) != bool:
            raise ValueError("The given parameter is not a boolean.")
        self.exclusiveMax = b

        msg = ""
        if self.__min > -1.e99 and self.__max == 1.e99:
            msg = self.tr("The float value must be greater than or equal to %.1f" % self.__min)
        elif self.__min == -1.e99 and self.__max < 1.e99:
            msg = self.tr("The float value must be lower than %.1f" % self.__max)
        elif self.__min > -1.e99 and self.__max < 1.e99:
            msg = self.tr("The float value must be greater than or equal to %.1f and lower than %.1f" % (self.__min, self.__max))

        self.parent.setStatusTip(str(msg))


    def setFixup(self, v):
        if type(v) != float:
            raise ValueError("The given parameter is not a float.")
        self.default = v
        self.fix = True


    def fixup(self, stri):
        if self.fix:
            if not stri:
                stri = str(self.default)


    def validate(self, stri, pos):
        """
        Validation method.

        QValidator.Invalid       0  The string is clearly invalid.
        QValidator.Intermediate  1  The string is a plausible intermediate value during editing.
        QValidator.Acceptable    2  The string is acceptable as a final result; i.e. it is valid.
        """
        state = QDoubleValidator.validate(self, stri, pos)[0]

        if state == QValidator.Acceptable:
            try:
                x = from_qvariant(stri, float)
            except Exception: # may be type error or localization issue
                x = 0.0
                state = QValidator.Intermediate

        if state == QValidator.Acceptable:
            if self.exclusiveMin and x == self.bottom():
                state = QValidator.Intermediate
            elif self.exclusiveMax and x == self.top():
                state = QValidator.Intermediate

        palette = self.parent.palette()

        if state != QValidator.Acceptable:
            palette.setColor(QPalette.Text, QColor("red"))
            self.parent.setPalette(palette)
        else:
            palette.setColor(QPalette.Text, QColor("black"))
            self.parent.setPalette(palette)

        self.state = state

        return (state, stri, pos)


    def tr(self, text):
        """
        """
        return text


class RegExpValidator(QRegExpValidator):
    """
    Validator for regular expression.
    """
    def __init__(self, parent, rx, forbidden_labels=None):
        """
        Initialization for validator
        """
        QRegExpValidator.__init__(self, parent)
        self.parent = parent
        self.state = QRegExpValidator.Invalid
        self.forbidden = forbidden_labels

        self.__validator = QRegExpValidator(rx, parent)

        if "{1," + str(LABEL_LENGTH_MAX) + "}" in rx.pattern():
            msg = self.tr("The maximum length of the label is %i characters" % LABEL_LENGTH_MAX)
            self.parent.setStatusTip(str(msg))


    def validate(self, stri, pos):
        """
        Validation method.

        QValidator.Invalid       0  The string is clearly invalid.
        QValidator.Intermediate  1  The string is a plausible intermediate value during editing.
        QValidator.Acceptable    2  The string is acceptable as a final result; i.e. it is valid.
        """
        state = self.__validator.validate(stri, pos)[0]

        if self.forbidden:
            if stri in self.forbidden:
                state = QValidator.Intermediate

        palette = self.parent.palette()

        if state == QValidator.Intermediate:
            palette.setColor(QPalette.Text, QColor("red"))
            self.parent.setPalette(palette)
        else:
            palette.setColor(QPalette.Text, QColor("black"))
            self.parent.setPalette(palette)

        self.state = state

        return (state, stri, pos)


    def tr(self, text):
        """
        """
        return text

#-------------------------------------------------------------------------------
# SpinBox progressing by multiplication and division
#-------------------------------------------------------------------------------

class RankSpinBoxWidget(QSpinBox):
    """
    Special Spin box for rank stepping.
    """
    def __init__(self, parent):
        """
        Constructor
        """
        QSpinBox.__init__(self, parent)

    def stepBy(self, steps):
        v = self.value()
        if steps > 0:
            self.setValue(v*2)
        elif steps < 0 and v > 1:
            self.setValue(v/2)

    def stepEnabled(self):
        v = self.value()
        if v < 2:
            return QAbstractSpinBox.StepUpEnabled
        else:
            return QAbstractSpinBox.StepUpEnabled | QAbstractSpinBox.StepDownEnabled

#-------------------------------------------------------------------------------
# SpinBox progressing by multiplication and division for Buffer size
#-------------------------------------------------------------------------------

class BufferSpinBoxWidget(QSpinBox):
    """
    Special Spin box for buffer size.
    """
    def __init__(self, parent):
        """
        Constructor
        """
        QSpinBox.__init__(self, parent)
        self.basesize = 1024*1024

    def stepBy(self, steps):
        v = self.value()
        if steps > 0:
            if v > 0:
                self.setValue(v*2)
            else:
                self.setValue(self.basesize)
        elif steps < 0 and v > 0:
            self.setValue(v/2)

    def textFromValue(self, v):
        """
        Define text to be shown.
        This text uses a local suffix (not that of the QSpinBox),
        as the suffix and value shown are dynamically related.
        """
        tv = v
        suffix = ''
        if v >= 1073741824 and v % 1073741824 == 0:
            tv = v // 1073741824
            suffix = ' GiB'
        elif v >= 1048576 and v % 1048576 == 0:
            tv = v // 1048576
            suffix = ' MiB'
        elif v >= 1024 and v % 1024 == 0:
            tv = v // 1024
            suffix = ' KiB'
        elif v > 0:
            tv = v
            suffix = ' B'
        else:
            tv = 0
            suffix = ''
        return QSpinBox.textFromValue(self, tv) + suffix

    def stepEnabled(self):
        v = self.value()
        if v < 1:
            return QAbstractSpinBox.StepUpEnabled
        else:
            return QAbstractSpinBox.StepUpEnabled | QAbstractSpinBox.StepDownEnabled

#-------------------------------------------------------------------------------
# Delegates to use
#-------------------------------------------------------------------------------

class LabelDelegate(QItemDelegate):
    """
    Delegate for lines
    """

    def __init__(self, parent=None, xml_model=None,
                 forbidden_labels=None,
                 accepted_regex=None,
                 auto_completion=[]):
        super(LabelDelegate, self).__init__(parent)

        self.parent = parent
        self.mdl    = xml_model

        self._forbidden_labels = forbidden_labels

        # Regex
        rx = accepted_regex
        if rx is None:
            rx = "[ -~]*"
#            "[\-_A-Za-z0-9]{1," + str(LABEL_LENGTH_MAX) + "}"
        self.regExp = QRegExp(rx)

        # Auto completion
        self._comp_list = auto_completion

    def createEditor(self, parent, option, index):

        editor = QLineEdit(parent)

        v = RegExpValidator(editor, self.regExp, self._forbidden_labels)
        editor.setValidator(v)

        # Autocompletion if necessary:
        if len(self._comp_list) > 0:
            completer = QCompleter()
            editor.setCompleter(completer)
            mc = QStringListModel()
            completer.setModel(mc)
            mc.setStringList(self._comp_list)

        # return editor
        return editor

    def setEditorData(self, editor, index):

        editor.setAutoFillBackground(True)

        v = from_qvariant(index.model().data(index, Qt.DisplayRole),
                          to_text_string)
        self.p_value = str(v)

        editor.setText(v)

    def setModelData(self, editor, model, index):

        if not editor.isModified():
            return

        if editor.validator().state == QValidator.Acceptable:
            p_value = str(editor.text())
            model.setData(index, p_value, Qt.DisplayRole)

class FloatDelegate(QItemDelegate):

    def __init__(self, parent=None, xml_model=None,
                 minVal=-1.e99, maxVal=+1.e99):

        super(FloatDelegate, self).__init__(parent)

        self.parent = parent
        self.mdl    = xml_model

        self._min   = minVal
        if type(self._min) != float:
            self._min = -1.e99

        self._max   = maxVal
        if type(self._max) != float:
            self._max = +1.e99

    def createEditor(self, parent, option, index):

        editor = QLineEdit(parent)

        validator = DoubleValidator(editor,
                                    min=self._min,
                                    max=self._max)

        editor.setValidator(validator)

        return editor

    def setEditorData(self, editor, index):

        editor.setAutoFillBackground(True)

        value = from_qvariant(index.model().data(index, Qt.DisplayRole),
                              to_text_string)
        editor.setText(value)

    def setModelData(self, editor, model, index):

        if editor.validator().state == QValidator.Acceptable:
            value = from_qvariant(editor.text(), float)
            selectionModel = self.parent.selectionModel()

            for idx in selectionModel.selectedIndexes():
                if idx.column() == index.column():
                    model.setData(idx, value, Qt.DisplayRole)

class IntegerDelegate(QItemDelegate):

    def __init__(self, parent=None, xml_model=None,
                 minVal=-10000000000, maxVal=+10000000000):

        super(IntegerDelegate, self).__init__(parent)

        self.parent  = parent
        self.mdl     = xml_model

        self._min = minVal
        if type(self._min) != int:
            self._min = -10000000000
        self._max = maxVal
        if type(self._max) != int:
            self._max = +10000000000


    def createEditor(self, parent, option, index):

        editor = QLineEdit(parent)

        validator = IntValidator(editor,
                                 min=self._min,
                                 max=self._max)

        editor.setValidator(validator)

        return editor

    def setEditorData(self, editor, index):

        editor.setAutoFillBackground(True)

        value = from_qvariant(index.model().data(index, Qt.DisplayRole),
                              to_text_string)
        editor.setText(value)

    def setModelData(self, editor, model, index):

        value = from_qvariant(editor.text(), int)

        if editor.validator().state == QValidator.Acceptable:
            selectionModel = self.parent.selectionModel()
            for idx in selectionModel.selectedIndexes():
                if idx.column() == index.column():
                    model.setData(idx, value, Qt.DisplayRole)


class ComboDelegate(QItemDelegate):

    def __init__(self, parent=None, xml_model=None,
                 opts_list=[], opts_state=None):

        super(ComboDelegate, self).__init__(parent)

        self.parent       = parent
        self.mdl          = xml_model

        if opts_state:
            if type(opts_state) != list:
                raise Exception("Wrong type for opts_state")
            if len(opts_state) != len(opts_list):
                raise Exception("Wrong length of opts_state")
        else:
            # To simplify the code which follows, we ensure that
            # opts_state is a list with the correct length
            opts_state = ["on"]*len(opts_list)

        self.opts_list = []
        for i, opt in enumerate(opts_list):

            if opts_state[i] not in ["on", "na", "off"]:
                msg="Wrong state for opts %s : %s" % (opt,opts_state[i])
                raise Exception(msg)

            ac = True
            av = True
            if opts_state[i] == "na":
                av = False
            elif opts_state[i] == "off":
                ac = False
                av = False

            self.opts_list.append({'name':opt,
                                   'active':ac,
                                   'available':av})


    def createEditor(self, parent, option, index):

        editor = QComboBox(parent)
        for opt in self.opts_list:
            name     = opt['name']
            isActive = opt['active']
            isAvail  = opt['available']

            editor.addItem(name)
            idx = editor.findText(name)
            editor.model().item(idx).setEnabled(isActive)
            if not isAvail:
                editor.setItemData(idx, QColor(Qt.red), Qt.TextColorRole)

        editor.installEventFilter(self)

        return editor


    def setEditorData(self, comboBox, index):

        string = from_qvariant(index.model().data(index, Qt.DisplayRole),
                               to_text_string)
        comboBox.setEditText(string)

    def setModelData(self, comboBox, model, index):

        value = comboBox.currentText()

        selectionModel = self.parent.selectionModel()

        for idx in selectionModel.selectedIndexes():
            if idx.column() == index.column():
                model.setData(idx, value, Qt.DisplayRole)


class BasicTableModel(QAbstractTableModel):

    def __init__(self, parent=QModelIndex(), xml_model=None, data=[[]], headers=[], default_row=[]):
        super(BasicTableModel, self).__init__(parent)
        self.parent = parent
        self.xml_model = xml_model
        self.data_table = data
        self.headers = headers
        self.default_row = default_row

    def rowCount(self, parent=QModelIndex()):
        return len(self.data_table)

    def columnCount(self, parent=QModelIndex()):
        return len(self.headers)

    def data(self, index, role=Qt.DisplayRole):
        if not index.isValid():
            return QVariant()
        elif role != Qt.DisplayRole:
            return QVariant()
        return self.data_table[index.row()][index.column()]

    def headerData(self, section, orientation, role=Qt.DisplayRole):
        if role == Qt.DisplayRole and orientation == Qt.Horizontal:
            return self.headers[section]
        return None

    def flags(self, index):
        if not index.isValid():
            return Qt.NoItemFlags
        else:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable

    def setData(self, index, value, role=Qt.EditRole):
        if (role == Qt.EditRole) and (index.isValid()):
            self.data_table[index.row()][index.column()] = value
            self.dataChanged.emit(index, index)
            return True
        return QAbstractTableModel.setData(index, value, role)

    def insertRows(self, position, rows=1, index=QModelIndex()):
        self.beginInsertRows(QModelIndex(), position, position + rows - 1)
        for row in range(rows):
            self.data_table.insert(position + row, self.default_row)
        self.endInsertRows()
        return True

    def removeRows(self, position, rows=1, index=QModelIndex()):
        self.beginRemoveRows(QModelIndex(), position, position + rows - 1)
        del self.data_table[position:position + rows]
        self.endRemoveRows()
        return True

# -------------------------------------------------------------------------------
# End of QtPage
# -------------------------------------------------------------------------------
