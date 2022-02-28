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
- QFileEditor
"""

#-------------------------------------------------------------------------------
# Standard modules
#-------------------------------------------------------------------------------

import sys, os, shutil
from code_saturne.gui.base import QtGui, QtCore, QtWidgets

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

import code_saturne.gui.base.resource_base_rc

#-------------------------------------------------------------------------------
# Local constants
#-------------------------------------------------------------------------------

_tab_size = 2

#-------------------------------------------------------------------------------
# Local functions and/or definitions
#-------------------------------------------------------------------------------

def loc_format(color, style=''):
    """
    Returns a TextCharFormat with the proper attributes
    """

    c = QtGui.QColor()
    c.setNamedColor(color)

    f = QtGui.QTextCharFormat()
    f.setForeground(c)

    # Bold font
    if 'bold' in style:
        f.setFontWeight(QtGui.QFont.Bold)

    # Italic font
    if 'italic' in style:
        f.setFontItalic(True)

    return f

format_styles = {'keyword'    : loc_format('blue', 'bold'),
                 'operator'   : loc_format('red', 'bold'),
                 'brace'      : loc_format('orange', 'bold'),
                 'string'     : loc_format('magenta', 'italic'),
                 'comment'    : loc_format('darkGreen', 'italic'),
                 'expression' : loc_format('black')}

#-------------------------------------------------------------------------------
# HighlightingRule class
#-------------------------------------------------------------------------------

class HighlightingRule():

    # ---------------------------------------------------------------
    def __init__(self, pattern, format):

        self.pattern = pattern
        self.format  = format
    # ---------------------------------------------------------------


#-------------------------------------------------------------------------------
# CodeEditor with line numbering
#-------------------------------------------------------------------------------

class LineNumberArea(QtWidgets.QWidget):

    def __init__(self, editor):
        # Handle the python2/python3 differences for super
        super(LineNumberArea, self).__init__(editor)

        self.editor = editor

    def sizeHint(self):
        return QtCore.QSize(self.editor.lineNumberAreaWidth(),0)

    def paintEvent(self, event):
        self.editor.lineNumberAreaPaintEvent(event)

class CodeEditor(QPlainTextEdit):
    def __init__(self):
        # Handle the python2/python3 differences for super
        super(CodeEditor, self).__init__()

        self.lineNumberArea = LineNumberArea(self)

        self.blockCountChanged.connect(self.updateLineNumberAreaWidth)
        self.updateRequest.connect(self.updateLineNumberArea)
        self.cursorPositionChanged.connect(self.highlightCurrentLine)

        self.updateLineNumberAreaWidth(0)


    def lineNumberAreaWidth(self):
        digits = 1
        count = max(1, self.blockCount())
        while count >= 10:
            count /= 10
            digits += 1
        space = 3 + self.fontMetrics().width('9') * digits
        return space


    def updateLineNumberAreaWidth(self, _):
        self.setViewportMargins(self.lineNumberAreaWidth(), 0, 0, 0)


    def updateLineNumberArea(self, rect, dy):

        if dy:
            self.lineNumberArea.scroll(0, dy)
        else:
            self.lineNumberArea.update(0, rect.y(), self.lineNumberArea.width(),
                       rect.height())

        if rect.contains(self.viewport().rect()):
            self.updateLineNumberAreaWidth(0)


    def resizeEvent(self, event):
        # Handle the python2/python3 differences for super
        super(CodeEditor, self).resizeEvent(event)

        cr = self.contentsRect();
        self.lineNumberArea.setGeometry(QtCore.QRect(cr.left(), cr.top(),
                    self.lineNumberAreaWidth(), cr.height()))


    def lineNumberAreaPaintEvent(self, event):
        mypainter = QtGui.QPainter(self.lineNumberArea)

        mypainter.fillRect(event.rect(), QtCore.Qt.lightGray)

        block = self.firstVisibleBlock()
        blockNumber = block.blockNumber()
        top = self.blockBoundingGeometry(block).translated(self.contentOffset()).top()
        bottom = top + self.blockBoundingRect(block).height()

        # Just to make sure I use the right font
        height = self.fontMetrics().height()
        while block.isValid() and (top <= event.rect().bottom()):
            if block.isVisible() and (bottom >= event.rect().top()):
                number = str(blockNumber + 1)
                mypainter.setPen(QtCore.Qt.black)
                mypainter.drawText(0, top, self.lineNumberArea.width(), height,
                 QtCore.Qt.AlignRight, number)

            block = block.next()
            top = bottom
            bottom = top + self.blockBoundingRect(block).height()
            blockNumber += 1


    def highlightCurrentLine(self):
        extraSelections = []

        if not self.isReadOnly():
            selection = QTextEdit.ExtraSelection()

            lineColor = QtGui.QColor(QtCore.Qt.yellow).lighter(160)

            selection.format.setBackground(lineColor)
            selection.format.setProperty(QtGui.QTextFormat.FullWidthSelection, True)
            selection.cursor = self.textCursor()
            selection.cursor.clearSelection()
            extraSelections.append(selection)
        self.setExtraSelections(extraSelections)

#-------------------------------------------------------------------------------
# QtextHighlighter class
#-------------------------------------------------------------------------------

class QtextHighlighter(QtGui.QSyntaxHighlighter):
    """
    Syntax highighting
    """

    # ---------------------------------------------------------------
    def __init__(self, parent, extension):

        QtGui.QSyntaxHighlighter.__init__(self, parent)
        self.parent = parent
        self.highlightingRules = []

        # Keywords (C or Fortran)
        fortran_kw = ['if', 'else', 'endif', 'do', 'enddo', 'end',
                      'implicit none', 'use', 'subroutine', 'function',
                      'double precision', 'real', 'integer', 'char',
                      'allocatable', 'allocate', 'deallocate', 'dimension',
                      'select case', 'call']

        c_kw       = ['if', 'else', 'for', 'switch', 'while',
                      '\#', 'include', 'pass', 'return', 'del', 'delete',
                      'assert', 'true', 'false', 'continue', 'break',
                      'fprintf', 'bft_printf', 'bft_printf_flush', 'bft_error',
                      'cs_real_t', 'cs_lnum_t', 'cs_real_3_t', 'int', 'char',
                      'string', 'void', 'double', 'const',
                      'BEGIN_C_DECLS', 'END_C_DECLS']

        py_kw      = ['if', 'elif', 'for', 'range', 'while', 'return', 'def',
                      'True', 'False']

        self.kw = []
        # Fortran
        if extension in ['f90', 'F90', 'F', 'f77']:
            for kw in fortran_kw:
                self.kw.append(kw)
                self.kw.append(kw.upper())
        # C/C++
        elif extension in ['c', 'cpp', 'cxx', 'c++']:
            for kw in c_kw:
                self.kw.append(kw)
                self.kw.append(kw.upper())
        # Python
        elif extension == 'py':
            for kw in py_kw:
                self.kw.append(kw)


        # Operators
        self.op = ['=', '==', '!=', '<', '>', '<=', '>=',
                   '\+', '-', '\*', '/', '\%', '\*\*',
                   '\+=', '-=', '\*=', '/=', '->', '=>',
                   '\^', '\|', '\&', '\|\|', '\&\&']

        # Braces
        self.br = ['\(', '\)', '\{', '\}', '\[', '\]']

        # RULES
        for kw in self.kw:
            p    = QtCore.QRegExp("\\b"+kw+ '\\b')
            rule = HighlightingRule(p, format_styles['keyword'])
            self.highlightingRules.append(rule)

        for op in self.op:
            p    = QtCore.QRegExp(op)
            rule = HighlightingRule(p, format_styles['operator'])
            self.highlightingRules.append(rule)

        for br in self.br:
            p    = QtCore.QRegExp(br)
            rule = HighlightingRule(p, format_styles['brace'])
            self.highlightingRules.append(rule)

        # strings
        ps = QtCore.QRegExp('"[^"\\]*(\\.[^"\\]*)*"')
        rs = HighlightingRule(ps, format_styles['string'])
        self.highlightingRules.append(rs)

        # comments
        pc = QtCore.QRegExp('//[^\n]*')
        rc = HighlightingRule(pc, format_styles['comment'])
        self.highlightingRules.append(rc)

        pcf = QtCore.QRegExp('![^\n]*')
        rcf = HighlightingRule(pcf, format_styles['comment'])
        self.highlightingRules.append(rcf)

        # numerals
        pn1 = QtCore.QRegExp('[+-]?[0-9]+[lL]?')
        rn1 = HighlightingRule(pn1, format_styles['expression'])
        self.highlightingRules.append(rn1)
        pn2 = QtCore.QRegExp('[+-]?0[xX][0-9A-Fa-f]+[lL]?')
        rn2 = HighlightingRule(pn2, format_styles['expression'])
        self.highlightingRules.append(rn2)
        pn3 = QtCore.QRegExp('[+-]?[0-9]+(?:\.[0-9]+)?(?:[eE][+-]?[0-9]+)?')
        rn3 = HighlightingRule(pn3, format_styles['expression'])
        self.highlightingRules.append(rn3)
    # ---------------------------------------------------------------


    # ---------------------------------------------------------------
    def highlightBlock(self, text):
        """
        Apply the syntax highlighting
        """
        for rule in self.highlightingRules:
            exp   = QtCore.QRegExp(rule.pattern)
            index = exp.indexIn(text)

            while index >= 0:
                length = exp.matchedLength()
                ok_to_highlight = True
                if len(text) > index+length:
                    if text[index+length] not in self.op+[' ']:
                        ok_to_highlight = False
                if text[index:index+length] not in self.op+self.br:
                    ok_to_highlight = True

                if ok_to_highlight:
                    self.setFormat(index, length, rule.format)
                if has_qstring:
                    index = text.indexOf(exp, index + length)
                else:
                    index = text.find(exp.cap(), index + length)

        self.setCurrentBlockState(0)

        # C/C++ comments
        self.highlightCommentsOverLines(text, "/\\*", "\\*/")
    # ---------------------------------------------------------------


    # ---------------------------------------------------------------
    def highlightCommentsOverLines(self, text, dls, dle):

        startExpression = QtCore.QRegExp(dls)
        endExpression   = QtCore.QRegExp(dle)
        ref_state = 1

        if self.previousBlockState() == ref_state:
            start = 0
            add   = 0

        else:
            start = startExpression.indexIn(text)
            add   = startExpression.matchedLength()


        while start >= 0:
            end = endExpression.indexIn(text, start + add)

            if end >= add:
                length = end - start + add + endExpression.matchedLength()
                self.setCurrentBlockState(0)

            else:
                self.setCurrentBlockState(ref_state)
                if has_qstring:
                    length = text.length() - start + add
                else:
                    length = len(text) - start + add

            self.setFormat(start, length, format_styles['comment'])
            start = endExpression.indexIn(text, start + length)
    # ---------------------------------------------------------------


#-------------------------------------------------------------------------------
# QMessageBox which expands
#-------------------------------------------------------------------------------

class QExpandingMessageBox(QMessageBox):
    """
    A QMessageBox which expands.
    """

    def __init__(self, parent=None):
        QMessageBox.__init__(self,parent=parent)
        self.setSizeGripEnabled(True)

    def event(self, ev):

        result = QMessageBox.event(self, ev)

        self.setMinimumHeight(10)
        self.setMaximumHeight(16777215)
        self.setMinimumWidth(10)
        self.setMaximumWidth(16777215)
        self.setSizePolicy(QSizePolicy.Expanding,
                           QSizePolicy.Expanding)

        text = self.findChild(QTextEdit)
        if text != None:
            self.setMinimumHeight(10)
            self.setMaximumHeight(16777215)
            self.setMinimumWidth(1050)
            self.setMaximumWidth(16777215)

            text.setMinimumHeight(10)
            text.setMaximumHeight(16777215)
            text.setMinimumWidth(1000)
            text.setMaximumWidth(16777215)
            text.setSizePolicy(QSizePolicy.Expanding,
                               QSizePolicy.Expanding)

        return result

#-------------------------------------------------------------------------------
# QFileEditor class
#-------------------------------------------------------------------------------

class FormWidget(QtWidgets.QWidget):
    """
    Main widget used to include both the browser and the editor zone
    """

    # ---------------------------------------------------------------
    def __init__(self, parent, wlist):
        super(FormWidget, self).__init__(parent)

        self.layout = QtWidgets.QGridLayout(self)

        n = len(wlist) - 1
        for i, w in enumerate(wlist):
            if i < n:
                w.setMaximumWidth(400)
                self.layout.addWidget(w, i, 0)
            else:
                self.layout.addWidget(w, 0, 1, 2, 1)

        self.setLayout(self.layout)
    # ---------------------------------------------------------------


#-------------------------------------------------------------------------------
# QFileSystemModel with modified header
#-------------------------------------------------------------------------------

class FileSystemModel(QtWidgets.QFileSystemModel):

    def __init__(self, title):
        """
        """
        QtWidgets.QFileSystemModel.__init__(self)
        self.title = title


    def headerData(self, section, orientation, role):
        if orientation == QtCore.Qt.Horizontal and role == QtCore.Qt.DisplayRole:
            if section == 0:
                return self.tr(self.title)
        return None


#-------------------------------------------------------------------------------
# Explorer class
#-------------------------------------------------------------------------------

class Explorer():
    """
    Editor class. Used for file editing and/or viewing
    """

    # ---------------------------------------------------------------
    def __init__(self, parent=None, root_dir=None, dir_type=None,
                 case_name=None, readOnly=False):

        self.parent = parent

        self.root_dir = root_dir
        self.dir_type = dir_type

        self.readOnly = readOnly
        self.readerMode = readOnly

        # Explorer
        self.explorer = self._initFileExplorer()
        self._initExplorerActions(case_name)

    # ---------------------------------------------------------------


    # ---------------------------------------------------------------
    def _initFileExplorer(self):
        """
        Create the File explorer object based on the QFileSystemModel widget.
        """

        if self.dir_type == 'SHARE':
            name = 'Reference'
        elif self.dir_type in ('SRC', 'DATA'):
            name = 'User files'
        else:
            name = 'Name'

        model = FileSystemModel(name)
        if self.root_dir:
            model.setRootPath(self.root_dir)

        tree = QtWidgets.QTreeView(None)

        tree.setModel(model)
        tree.setSortingEnabled(True)
        tree.setWindowTitle('Explorer')
        if self.root_dir:
            tree.setRootIndex(model.index(self.root_dir))

        # Hide unnecessary columns
        nc = tree.header().count()

        for i in range(1, nc):
            tree.hideColumn(i)

        # Right click menu
        tree.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        tree.customContextMenuRequested.connect(self.explorerContextMenu)

        # Double click
        tree.doubleClicked.connect(self._explorerDoubleClick)

        return tree;
    # ---------------------------------------------------------------


    # ---------------------------------------------------------------
    def _initExplorerActions(self, case_name=None):
        """
        Create explorer actions dictionary
        """

        if case_name:
            case_dir_name = str(case_name)
        else:
            case_dir_name = 'SRC'

        _editAction = QAction(self.explorer.model())
        _editAction.setText('Edit file')
        _editAction.triggered.connect(self.parent._editSelectedFile)

        _viewAction = QAction(self.explorer.model())
        _viewAction.setText('View file')
        _viewAction.triggered.connect(self.parent._viewSelectedFile)

        _copyAction = QAction(self.explorer.model())
        _copyAction.setText('Copy to ' + case_dir_name)
        _copyAction.triggered.connect(self.parent._copySelectedFile)

        _removeAction = QAction(self.explorer.model())
        _removeAction.setText('Remove from ' + case_dir_name)
        _removeAction.triggered.connect(self.parent._removeSelectedFile)

        _restoreAction = QAction(self.explorer.model())
        _restoreAction.setText('Move to ' + case_dir_name)
        _restoreAction.triggered.connect(self.parent._restoreSelectedFile)

        _deleteAction = QAction(self.explorer.model())
        _deleteAction.setText('Delete')
        _deleteAction.triggered.connect(self.parent._deleteSelectedFile)

        self._explorerActions = {'edit': _editAction,
                                 'view': _viewAction,
                                 'copy': _copyAction,
                                 'remove': _removeAction,
                                 'restore': _restoreAction,
                                 'delete': _deleteAction}
    # ---------------------------------------------------------------


    # ---------------------------------------------------------------
    def _updateCurrentSelection(self):
        """
        Update the current selection
        """
        # Find file position (SRC, REFERENCE, EXAMPLES, other)
        path2file = ''
        for idx in self.explorer.selectedIndexes():
            fname = idx.data(QtCore.Qt.DisplayRole)
            c = idx
            p = c.parent()
            ps = p.data(QtCore.Qt.DisplayRole)
            while True:
                ctxt = c.data(QtCore.Qt.DisplayRole)
                ptxt = p.data(QtCore.Qt.DisplayRole)
                if ptxt in [None, self.parent.case_name]:
                    pe = ptxt
                    break
                path2file = ptxt + '/' + path2file
                c = p
                p = c.parent()

        self.parent._currentSelection = {'filename':fname,
                                         'subpath' :path2file,
                                         'filedir' :ps,
                                         'origdir' :pe}

        return
    # ---------------------------------------------------------------


    # ---------------------------------------------------------------
    def _explorerDoubleClick(self):
        """
        Double click action
        """

        self._updateCurrentSelection()

        clicked = os.path.join(self.parent._currentSelection['subpath'],
                               self.parent._currentSelection['filename'])

        # To ensure that os.path.isdir works correctly we use the full path
        # to the object which is selected in the menu
        if self.root_dir:
            clicked = os.path.join(self.root_dir, clicked)

        edit_list = ['SRC']

        if not os.path.isdir(clicked):
            if self.parent._currentSelection['filedir'] in edit_list:
                self.parent._editSelectedFile()
            else:
                self.parent._viewSelectedFile()

    # ---------------------------------------------------------------


    # ---------------------------------------------------------------
    def explorerContextMenu(self, position):
        """
        Custom menu for the mouse right-click.
        Depends on whether the file is in the SRC, SRC/subfolder
        or RESU/subfolder.
        Possible actions are 'edit', 'view' and 'copy' (to SRC)
        """

        self._updateCurrentSelection()

        path2file = self.parent._currentSelection['subpath']
        fname     = self.parent._currentSelection['filename']
        pe        = self.parent._currentSelection['origdir']
        ps        = self.parent._currentSelection['filedir']

        self._contextMenu = QMenu()

        if (path2file == '' or path2file is None ) and self.root_dir:
            path2file = self.root_dir

        if self.dir_type == 'SHARE':
            if not os.path.isdir(os.path.join(path2file, fname)):
                self._contextMenu.addAction(self._explorerActions['view'])
                self._contextMenu.addAction(self._explorerActions['copy'])
        elif pe == 'SRC':
            if not os.path.isdir(os.path.join(path2file, fname)):
                if ps == 'SRC':
                    self._contextMenu.addAction(self._explorerActions['edit'])
                    self._contextMenu.addAction(self._explorerActions['remove'])
                elif ps in ['EXAMPLES', 'REFERENCE']:
                    self._contextMenu.addAction(self._explorerActions['view'])
                    self._contextMenu.addAction(self._explorerActions['copy'])
                elif ps in ['DRAFT']:
                    self._contextMenu.addAction(self._explorerActions['view'])
                    self._contextMenu.addAction(self._explorerActions['restore'])
                    self._contextMenu.addAction(self._explorerActions['delete'])
        elif pe == 'DATA':
            if not os.path.isdir(os.path.join(path2file, fname)):
                if ps == 'DATA':
                    if fname not in ('setup.xml', 'run.cfg'):
                        self._contextMenu.addAction(self._explorerActions['edit'])
                        self._contextMenu.addAction(self._explorerActions['remove'])
                    else:
                        self._contextMenu.addAction(self._explorerActions['view'])
                elif ps in ['REFERENCE']:
                    self._contextMenu.addAction(self._explorerActions['view'])
                    self._contextMenu.addAction(self._explorerActions['copy'])
                elif ps in ['DRAFT']:
                    self._contextMenu.addAction(self._explorerActions['view'])
                    self._contextMenu.addAction(self._explorerActions['restore'])
                    self._contextMenu.addAction(self._explorerActions['delete'])
        else:
            if not os.path.isdir(os.path.join(path2file, fname)):
                self._contextMenu.addAction(self._explorerActions['view'])

        self._contextMenu.exec_(self.explorer.viewport().mapToGlobal(position))
    # ---------------------------------------------------------------


#-------------------------------------------------------------------------------
# QFileEditor class
#-------------------------------------------------------------------------------

class QFileEditor(QMainWindow):
    """
    Editor class. Used for file editing and/or viewing
    """

    # ---------------------------------------------------------------
    def __init__(self, parent=None, case_dir=None, reference_dir=None,
                 readOnly=False, noOpen=False, useHighlight=True):
        super(QFileEditor, self).__init__(parent)
        self.setGeometry(50, 50, 500, 300)

        self.setWindowTitle("code_saturne built-in file editor")
        self.parent = parent

        self.case_dir = case_dir
        if self.case_dir:
            self.case_name = os.path.split(case_dir)[-1]

        self.last_dir = case_dir

        self.readOnly = readOnly
        self.readerMode = readOnly

        # Activate text highlight
        self.useHighlight = useHighlight

        self.opened = False
        self.saved  = True

        # Open file action
        open_img_path = ":/icons/22x22/document-open.png"
        icon_open     = QtGui.QIcon()
        icon_open.addPixmap(QtGui.QPixmap(_fromUtf8(open_img_path)),
                            QtGui.QIcon.Normal,
                            QtGui.QIcon.Off)
        self.openFileAction = QAction(icon_open, "Open", self)
        self.openFileAction.setShortcut("Ctrl+O")
        self.openFileAction.setStatusTip('Open File')
        self.openFileAction.triggered.connect(self.openFileForAction)

        # New file action
        new_img_path = ":/icons/22x22/document-new.png"
        icon_new     = QtGui.QIcon()
        icon_new.addPixmap(QtGui.QPixmap(_fromUtf8(new_img_path)),
                          QtGui.QIcon.Normal,
                          QtGui.QIcon.Off)
        self.newFileAction = QAction(icon_new, "New", self)
        self.newFileAction.setShortcut("Ctrl+E")
        self.newFileAction.setStatusTip('Create new file')
        self.newFileAction.triggered.connect(self.newFile)

        # Save action
        save_img_path = ":/icons/22x22/document-save.png"
        icon_save     = QtGui.QIcon()
        icon_save.addPixmap(QtGui.QPixmap(_fromUtf8(save_img_path)),
                          QtGui.QIcon.Normal,
                          QtGui.QIcon.Off)
        self.saveFileAction = QAction(icon_save, "Save", self)
        self.saveFileAction.setShortcut("Ctrl+S")
        self.saveFileAction.setStatusTip('Save file')
        self.saveFileAction.triggered.connect(self.saveFile)

        # Save as action
        saveas_img_path = ":/icons/22x22/document-save-as.png"
        icon_saveas     = QtGui.QIcon()
        icon_saveas.addPixmap(QtGui.QPixmap(_fromUtf8(saveas_img_path)),
                              QtGui.QIcon.Normal,
                              QtGui.QIcon.Off)
        self.saveFileAsAction = QAction(icon_saveas, "Save as", self)
        self.saveFileAsAction.setStatusTip('Save file as')
        self.saveFileAsAction.triggered.connect(self.saveFileAs)

        # Close file action
        close_img_path = ":/icons/22x22/process-stop.png"
        icon_close     = QtGui.QIcon()
        icon_close.addPixmap(QtGui.QPixmap(_fromUtf8(close_img_path)),
                             QtGui.QIcon.Normal,
                             QtGui.QIcon.Off)
        self.closeFileAction = QAction(icon_close, "Close file", self)
        self.closeFileAction.setShortcut("Ctrl+Q")
        self.closeFileAction.setStatusTip('Close opened file')
        self.closeFileAction.triggered.connect(self.closeOpenedFile)

        # Exit editor action
        quit_img_path = ":/icons/22x22/system-log-out.png"
        icon_quit     = QtGui.QIcon()
        icon_quit.addPixmap(QtGui.QPixmap(_fromUtf8(quit_img_path)),
                          QtGui.QIcon.Normal,
                          QtGui.QIcon.Off)
        self.quitAction = QAction(icon_quit, "Quit", self)
        self.quitAction.setStatusTip('Quit the editor')
        self.quitAction.triggered.connect(self.closeApplication)

        self.statusBar()

        # File toolbar
        self.toolbar = self.addToolBar("Options")

        self.toolbar.addAction(self.newFileAction)
        if not noOpen:
            self.toolbar.addAction(self.openFileAction)
        self.toolbar.addAction(self.saveFileAction)
        self.toolbar.addAction(self.saveFileAsAction)
        self.toolbar.addAction(self.closeFileAction)
        self.toolbar.addAction(self.quitAction)

        # File menu
        self.mainMenu = self.menuBar()

        self.fileMenu = self.mainMenu.addMenu('&File')
        self.fileMenu.addAction(self.newFileAction)
        if not noOpen:
            self.fileMenu.addAction(self.openFileAction)
        self.fileMenu.addAction(self.saveFileAction)
        self.fileMenu.addAction(self.saveFileAsAction)
        self.fileMenu.addAction(self.closeFileAction)
        self.fileMenu.addAction(self.quitAction)

        # Explorer
        self.explorer = Explorer(parent=self,
                                 root_dir=self.case_dir,
                                 dir_type=self.case_name,
                                 case_name=self.case_name)

        # Explorer
        self.explorer_ref = None
        if reference_dir:
            self.explorer_ref = Explorer(parent=self,
                                         root_dir=reference_dir,
                                         dir_type='SHARE',
                                         case_name=self.case_name)

        # Editor
        self.textEdit = self._initFileEditor()

        # Settings
        settings = QtCore.QSettings()

        try:
            # API 2
            self.restoreGeometry(settings.value("MainWindow/Geometry", QtCore.QByteArray()))
            self.restoreState(settings.value("MainWindow/State", QtCore.QByteArray()))
        except:
            # API 1
            self.recentFiles = settings.value("RecentFiles").toStringList()
            self.restoreGeometry(settings.value("MainWindow/Geometry").toByteArray())
            self.restoreState(settings.value("MainWindow/State").toByteArray())

        # file attributes
        self.filename = ""
        self.file_extension  = ""

        if self.explorer_ref:
            self.mainWidget = FormWidget(self, [self.explorer.explorer,
                                                self.explorer_ref.explorer,
                                                self.textEdit])
        else:
            self.mainWidget = FormWidget(self, [self.explorer.explorer,
                                                self.textEdit])

        self.setCentralWidget(self.mainWidget)
    # ---------------------------------------------------------------


    # ---------------------------------------------------------------
    def _initFileEditor(self):
        """
        Create the Editor widget based on QTextEdit
        """

        # Font
        base_font = QtGui.QFont()
        base_font.setFamily("Courier")
        base_font.setStyleHint(QtGui.QFont.Monospace)
        base_font.setFixedPitch(True)
        base_font.setPointSize(10)

        font_metrics = QtGui.QFontMetrics(base_font)
        _tab_string = ''
        for i in range(_tab_size):
            _tab_string += ' '

        # Main text zone
        textEdit = CodeEditor()
        textEdit.setFont(base_font)
        textEdit.textChanged.connect(self.updateFileState)
        textEdit.setReadOnly(self.readOnly)
        policy = textEdit.sizePolicy()
        policy.setHorizontalPolicy(QSizePolicy.Expanding)
        textEdit.setSizePolicy(policy)

        # tab
        textEdit.setTabStopWidth(font_metrics.width(_tab_string))

        return textEdit
    # ---------------------------------------------------------------


    # ---------------------------------------------------------------
    def _initFileExplorer(self, base_dir=None, name="User Files"):
        """
        Create the File explorer object based on the QFileSystemModel widget.
        """

        #model = QtWidgets.QFileSystemModel()
        model = FileSystemModel(name)
        rootp = ''
        if base_dir:
            rootp = base_dir
        elif self.case_dir:
            rootp = self.case_dir

        model.setRootPath(rootp)

        tree = QtWidgets.QTreeView(None)

        tree.setModel(model)
        tree.setSortingEnabled(True)
        tree.setWindowTitle('Explorer')
        tree.setRootIndex(model.index(rootp))

        # Hide unnecessary columns
        nc = tree.header().count()

        for i in range(1, nc):
            tree.hideColumn(i)

        # Right click menu
        tree.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        tree.customContextMenuRequested.connect(self.explorerContextMenu)

        # Double click
        tree.doubleClicked.connect(self._explorerDoubleClick)

        return tree;
    # ---------------------------------------------------------------


    # ---------------------------------------------------------------
    def _editSelectedFile(self):
        """
        Edit action for mouse right-click
        """

        self.readOnly = False

        t = "Editor: %s" % (self._currentSelection['filename'])
        self.setWindowTitle(t)

        fn = os.path.join(self.case_dir,
                          self._currentSelection['subpath'],
                          self._currentSelection['filename'])
        self.openFile(fn=fn)
    # ---------------------------------------------------------------


    # ---------------------------------------------------------------
    def _viewSelectedFile(self):
        """
        View action for mouse left-click
        """

        self.readOnly = True

        t = "Viewer: %s" % (self._currentSelection['filename'])
        self.setWindowTitle(t)

        fn = os.path.join(self.case_dir,
                          self._currentSelection['subpath'],
                          self._currentSelection['filename'])
        self.openFile(fn=fn)
    # ---------------------------------------------------------------


    # ---------------------------------------------------------------
    def _copySelectedFile(self):
        """
        Copy files in subdirectories, such as REFERENCES or EXAMPLES
        to the SRC folder. Used by the mouse right-click
        """

        src_path = os.path.join(self.case_dir,
                                self._currentSelection['subpath'],
                                self._currentSelection['filename'])

        if self.case_name in ('SRC', 'DATA'):
            trg_path = os.path.join(self.case_dir,
                                    self._currentSelection['filename'])
        else:
            sp = self._currentSelection['subpath']
            while '/' in sp and len(sp) > 3:
                e1, e2 = os.path.split(sp)
                if e2 in ('SRC', 'DATA'):
                    break
                else:
                    sp = e1

            trg_path = os.path.join(sp, self._currentSelection['filename'])

        shutil.copy2(src_path, trg_path)
    # ---------------------------------------------------------------


    # ---------------------------------------------------------------
    def _removeSelectedFile(self):
        """
        Remove a file from the SRC dir
        """

        title = "Remove file"
        question = "Remove %s from the SRC folder (Stored in DRAFT) ?" % (self._currentSelection['filename'])

        choice = QMessageBox.question(self,
                                      title,
                                      question,
                                      QMessageBox.Yes | QMessageBox.No)

        if choice == QMessageBox.Yes:
            fn = os.path.join(self.case_dir,
                              self._currentSelection['subpath'],
                              self._currentSelection['filename'])


            draft = os.path.join(self.case_dir,
                               self._currentSelection['subpath'],
                               'DRAFT')
            if not os.path.exists(draft):
                os.mkdir(draft)
            fn2 = os.path.join(draft, self._currentSelection['filename'])

            if os.path.exists(fn2):
                q = 'A file named %s allready exists in DRAFT.\nDo you want to overwrite it?' % (self._currentSelection['filename'])
                choice2 = QMessageBox.question(self,
                                                     '',
                                                     q,
                                                     QMessageBox.Yes | QMessageBox.No)
                if choice2 == QMessageBox.No:
                    return

            shutil.move(fn, fn2)
        else:
            pass
    # ---------------------------------------------------------------


    # ---------------------------------------------------------------
    def _restoreSelectedFile(self):
        """
        Move a file from DRAFT to the SRC folder
        """

        title = "Move to SRC"
        question = "Move file %s from DRAFT to SRC folder ?" % (self._currentSelection['filename'])

        choice = QMessageBox.question(self,
                                      title,
                                      question,
                                      QMessageBox.Yes | QMessageBox.No)

        if choice == QMessageBox.Yes:
            fn = os.path.join(self.case_dir,
                              self._currentSelection['subpath'],
                              self._currentSelection['filename'])

            fn2 = os.path.join(self.case_dir, self._currentSelection['filename'])

            if os.path.exists(fn2):
                q = 'A file named %s allready exists in SRC\nDo you want to overwrite it?' % (self._currentSelection['filename'])
                choice2 = QMessageBox.question(self, '', q,
                                               QMessageBox.Yes | QMessageBox.No)

                if choice2 == QMessageBox.No:
                    return

            shutil.move(fn, fn2)

        else:
            pass
    # ---------------------------------------------------------------


    # ---------------------------------------------------------------
    def _deleteSelectedFile(self):
        """
        Remove a file from the SRC dir
        """

        title = "Delete file"
        question = "Really delete %s ?" % (self._currentSelection['filename'])

        choice = QMessageBox.question(self,
                                      title,
                                      question,
                                      QMessageBox.Yes | QMessageBox.No)

        if choice == QMessageBox.Yes:
            fn = os.path.join(self.case_dir,
                              self._currentSelection['subpath'],
                              self._currentSelection['filename'])

            try:
                os.remove(fn)
            except Exception:
                # TODO add error popup
                pass

            d = os.path.split(fn)[0]
            if os.path.basename(d) in ('DRAFT', 'STASH'):
                l = os.listdir(d)
                if len(l) < 1:
                    try:
                        os.rmdir(d)
                    except Exception:
                        pass
        else:
            pass
    # ---------------------------------------------------------------


    # ---------------------------------------------------------------
    def updateFileState(self, new_state = False):
        """
        Update file state (saved or not)
        """
        self.saved  = new_state
        # To ensure syntax highlighting while modifying the text
        self.textEdit.viewport().update()
    # ---------------------------------------------------------------


    # ---------------------------------------------------------------
    def openFile(self, fn = None):
        """
        Open a file in the editor
        """

        if not self.saved:
            self.closeOpenedFile()

        if fn:
            self.filename = fn
        else:
            self.filename = QFileDialog.getOpenFileName(self, 'Open File',
                                                        self.last_dir)

        if self.filename:
            self.last_dir = os.path.split(self.filename)[0]

        self.textEdit.setReadOnly(self.readOnly)
        self.saveFileAction.setEnabled(not self.readOnly)

        if self.filename != None and self.filename != '':
            file = open(self.filename, 'r')
            self.file_extension = self.filename.split('.')[-1]

            self.newFile()
            with file:
                text = file.read()
                self.textEdit.setPlainText(text)
                self.updateFileState(True)
    # ---------------------------------------------------------------


    # ---------------------------------------------------------------
    def openFileForAction(self, fn = None):

        if self.readOnly != self.readerMode:
            self.readOnly = self.readerMode

        self.openFile(fn)
    # ---------------------------------------------------------------


    # ---------------------------------------------------------------
    def newFile(self):
        """
        Create a new file (blank)
        """

        self.opened = True
        self.updateFileState(False)
        if self.useHighlight:
            hl = QtextHighlighter(self.textEdit.document(), self.file_extension)
        self.textEdit.show()
    # ---------------------------------------------------------------


    # ---------------------------------------------------------------
    def saveFile(self):
        """
        Save file
        """
        if not self.opened:
            return

        if self.filename != None and self.filename != '':
            file = open(self.filename,'w')
            text = self.textEdit.toPlainText()
            file.write(text)
            file.close()

            self.updateFileState(True)

        else:
            self.saveFileAs()
    # ---------------------------------------------------------------


    # ---------------------------------------------------------------
    def saveFileAs(self):
        """
        Save file as
        """
        if not self.opened:
            return

        ret = QFileDialog.getSaveFileName(self, 'Save File')

        if type(ret) == str:
            self.filename = ret
        elif type(ret) == tuple:
            self.filename = ret[0]
        else:
            raise Exception("Uknown return type for 'QFileDialog.getSaveFileName'")

        self.last_dir = os.path.split(self.filename)[0]

        if self.filename != None and self.filename != '':
            file = open(self.filename,'w')
            text = self.textEdit.toPlainText()
            file.write(text)
            file.close()

            self.updateFileState(True)
    # ---------------------------------------------------------------


    # ---------------------------------------------------------------
    def closeOpenedFile(self):
        """
        Close an opened file
        """

        if self.saved == False and self.readOnly == False:
            choice = QMessageBox.question(self, 'Built-in editor',
                                          'File changed.\nDo you want to save?',
                                          QMessageBox.Yes | QMessageBox.No)
            if choice == QMessageBox.Yes:
                self.saveFile()
            else:
                pass

        self.saved  = True
        self.opened = False

        self.filename = ''
        self.textEdit.setPlainText('')
    # ---------------------------------------------------------------


    # ---------------------------------------------------------------
    def closeApplication(self):
        """
        Close the editor
        """
        if self.opened == True:
            choice = QMessageBox.question(self, 'Built-in editor',
                                          "Exit text editor?",
                                          QMessageBox.Yes | QMessageBox.No)
        else:
            choice = QMessageBox.Yes

        if choice == QMessageBox.Yes:
            self.closeOpenedFile()

            settings = QtCore.QSettings()
            settings.setValue("MainWindow/Geometry",
                              self.saveGeometry())

            self.close()
            return 0
        else:
            return 1
    # ---------------------------------------------------------------


    # ---------------------------------------------------------------
    def closeEvent(self, event):

        decision = self.closeApplication()
        if decision == 1:
            event.ignore()
    # ---------------------------------------------------------------


#-------------------------------------------------------------------------------
# END OF FILE
#-------------------------------------------------------------------------------
