# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2020 EDF S.A.
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

Py2 = sys.version[0] == '2'
Py3 = sys.version[0] == '3'


#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------

from code_saturne.Base.QtCore    import *
from code_saturne.Base.QtGui     import *
from code_saturne.Base.QtWidgets import *

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
if Py2:
    TEXT_TYPES = (str, unicode)
else:
    TEXT_TYPES = (str,)

#==============================================================================
# Strings
#==============================================================================

def is_text_string(obj):
    """Return True if `obj` is a text string,
              False if it is anything else,
                    like binary data (Python 3) or
                    QString (Python 2, PyQt API #1)"""
    if Py2:
        return isinstance(obj, basestring)
    else:
        return isinstance(obj, str)

def to_text_string(obj, encoding=None):
    """Convert `obj` to (unicode) text string"""
    if Py2:
        if encoding is None:
            return unicode(obj)
        else:
            return unicode(obj, encoding)
    else:
        if encoding is None:
            return str(obj)
        elif isinstance(obj, str):
            # In case this function is not used properly, this could happen
            return obj
        else:
            return str(obj, encoding)

#==============================================================================
# QVariant conversion utilities
#==============================================================================
PYQT_API_1 = False
import collections

if os.environ.get('QT_API', 'pyqt') == 'pyqt':
    import sip
    if QT_API == "PYQT4":
        try:
            PYQT_API_1 = sip.getapi('QVariant') == 1
        except AttributeError:
            PYQT_API_1 = True

    def to_qvariant(pyobj=None):
        """Convert Python object to QVariant
        This is a transitional function from PyQt API #1 (QVariant exist)
        to PyQt API #2 and Pyside (QVariant does not exist)"""
        if PYQT_API_1:
            from code_saturne.Base.QtCore import QVariant
            return QVariant(pyobj)
        else:
            return pyobj

    def from_qvariant(qobj=None, convfunc=None):
        """Convert QVariant object to Python object
        This is a transitional function from PyQt API #1 (QVariant exists)
        to PyQt API #2 and Pyside (QVariant does not exist)"""
        if PYQT_API_1:
            assert isinstance(convfunc, collections.Callable)
            if convfunc in TEXT_TYPES or convfunc is to_text_string:
                return convfunc(qobj.toString())
            elif convfunc is bool:
                return qobj.toBool()
            elif convfunc is int:
                return qobj.toInt()[0]
            elif convfunc is float:
                return qobj.toDouble()[0]
            else:
                return convfunc(qobj)
        else:
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
    def to_qvariant(obj=None):
        """Convert Python object to QVariant
        This is a transitional function from PyQt API#1 (QVariant exist)
        to PyQt API#2 and Pyside (QVariant does not exist)"""
        return obj

    def from_qvariant(qobj=None, pytype=None):
        """Convert QVariant object to Python object
        This is a transitional function from PyQt API #1 (QVariant exist)
        to PyQt API #2 and Pyside (QVariant does not exist)"""
        return qobj

def qbytearray_to_str(qba):
    """Convert QByteArray object to str in a way compatible with Python 2/3"""
    return str(bytes(qba.toHex().data()).decode())


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
        from code_saturne.Base.QtCore import QString
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

        self.items   = []
        self.combo.setModel(self.model)


    def addItem(self, str_view, str_model="", warn=False):
        """
        Insert an item in the model.

        str_view: string to be displayed in the view.
        For example, 'Eulerian/Lagrangian Multi-phase Treatment'

        str_model: correponding string used in the model.
        For example, 'lagrangian'

        warn: If True, entry is marked with a color.
        """
        item  = QStandardItem(str(str_view))

        index = self.last
        self.model.setItem(index, item)

        if warn:
            self.combo.setItemData(index,
                                   QColor(Qt.red),
                                   Qt.TextColorRole)

        self.last = index + 1

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
            self.dicoV2M[new_str_view]  = new_str_model


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
                print(str_model, " is not in list: ", str(self.items))
                print("Value reset to ", str(self.items[0]),
                      " but should be checked.")
                index = 0

        elif str_view:
            str_model = self.dicoV2M[str_view]
            index = self.items.index(str_model)
            self.combo.setCurrentIndex(index)


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

        if PYQT_API_1:
            return (state, pos)
        else:
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
            msg = self.tr("The float value must be greater than %.1f" % min)
        elif min == -1.e99 and max < 1.e99:
            msg = self.tr("The float value must be lower than %.1f" % max)
        elif min > -1.e99 and max < 1.e99:
            msg = self.tr("The float value must be between than %.1f and %.1f" % (min, max))

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

        if PYQT_API_1:
            return (state, pos)
        else:
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

        if PYQT_API_1:
            return (state, pos)
        else:
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
            tv = v / 1073741824
            suffix = ' GiB'
        elif v >= 1048576 and v % 1048576 == 0:
            tv = v / 1048576
            suffix = ' MiB'
        elif v >= 1024 and v % 1024 == 0:
            tv = v / 1024
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
# End of QtPage
#-------------------------------------------------------------------------------
