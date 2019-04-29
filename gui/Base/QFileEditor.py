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
- QFileEditor
"""

#-------------------------------------------------------------------------------
# Standard modules
#-------------------------------------------------------------------------------

import sys, os, shutil
from code_saturne.Base import QtGui, QtCore, QtWidgets

# Check if QString exists
has_qstring = True
try:
    from code_saturne.Base.QtCore import QString
    _fromUtf8 = QString.fromUtf8
except ImportError:
    has_qstring = False
    def _fromUtf8(s):
        return s

    def QString(s):
        return s

import resource_base_rc

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
# QtextHighlighter class
#-------------------------------------------------------------------------------

class QtextHighlighter(QtGui.QSyntaxHighlighter):
    """
    Syntax highighting
    """

    # ---------------------------------------------------------------
    def __init__(self, parent):

        QtGui.QSyntaxHighlighter.__init__(self, parent)
        self.parent = parent
        self.highlightingRules = []

        # Keywords (C or Fortran)
        self.kw = ['if', 'else', 'endif', '\#', 'include',
                   'void', 'int', 'integer', 'double', 'const',
                   'fprintf', 'bft_printf', 'bft_printf_flush',
                   'cs_real_t',
                   'subroutine', 'function', 'def',
                   'double precision', 'use', 'implicit none',
                   'allocatable', 'dimension', 'string', 'float',
                   'allocate', 'deallocate',
                   'char', 'for', 'while', 'assert',
                   'continue', 'break', 'switch',
                   'del', 'pass', 'return', 'true', 'false']

        # Operators
        self.op = ['=', '==', '!=', '<', '>', '<=', '>=',
                   '\+', '-', '\*', '/', '\%', '\*\*',
                   '\+=', '-=', '\*=', '/=',
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
# QFileEditor class
#-------------------------------------------------------------------------------
class FormWidget(QtWidgets.QWidget):
    """
    Main widget used to include both the browser and the editor zone
    """

    # ---------------------------------------------------------------
    def __init__(self, parent, wlist):
        super(FormWidget, self).__init__(parent)

        self.layout = QtWidgets.QHBoxLayout(self)

        for w in wlist:
            if w == wlist[0]:
                w.setMaximumWidth(400)
            self.layout.addWidget(w)

        self.setLayout(self.layout)
    # ---------------------------------------------------------------


#-------------------------------------------------------------------------------
# QFileEditor class
#-------------------------------------------------------------------------------

class QFileEditor(QtGui.QMainWindow):
    """
    Editor class. Used for file editing and/or viewing
    """

    # ---------------------------------------------------------------
    def __init__(self, parent=None, case_dir=None, readOnly=False):
        super(QFileEditor, self).__init__(parent)
        self.setGeometry(50, 50, 500, 300)

        self.setWindowTitle("Code_Saturne built-in file editor")
        self.parent = parent

        self.case_dir = case_dir
        if self.case_dir:
            self.case_name = os.path.split(case_dir)[-1]

        self.last_dir = case_dir

        self.readOnly = readOnly

        self.saved = True

        self.dialog = QtGui.QFileDialog(self)

        # Open file action
        open_img_path = ":/icons/22x22/document-open.png"
        icon_open     = QtGui.QIcon()
        icon_open.addPixmap(QtGui.QPixmap(_fromUtf8(open_img_path)),
                            QtGui.QIcon.Normal,
                            QtGui.QIcon.Off)
        self.openFileAction = QtGui.QAction(icon_open, "Open", self)
        self.openFileAction.setShortcut("Ctrl+O")
        self.openFileAction.setStatusTip('Open File')
        self.openFileAction.triggered.connect(self.openFile)

        # New file action
        new_img_path = ":/icons/22x22/document-new.png"
        icon_new     = QtGui.QIcon()
        icon_new.addPixmap(QtGui.QPixmap(_fromUtf8(new_img_path)),
                          QtGui.QIcon.Normal,
                          QtGui.QIcon.Off)
        self.newFileAction = QtGui.QAction(icon_new, "New", self)
        self.newFileAction.setShortcut("Ctrl+E")
        self.newFileAction.setStatusTip('Create new file')
        self.newFileAction.triggered.connect(self.newFile)

        # Save action
        save_img_path = ":/icons/22x22/document-save.png"
        icon_save     = QtGui.QIcon()
        icon_save.addPixmap(QtGui.QPixmap(_fromUtf8(save_img_path)),
                          QtGui.QIcon.Normal,
                          QtGui.QIcon.Off)
        self.saveFileAction = QtGui.QAction(icon_save, "Save", self)
        self.saveFileAction.setShortcut("Ctrl+S")
        self.saveFileAction.setStatusTip('Save file')
        self.saveFileAction.triggered.connect(self.saveFile)

        # Save as action
        saveas_img_path = ":/icons/22x22/document-save-as.png"
        icon_saveas     = QtGui.QIcon()
        icon_saveas.addPixmap(QtGui.QPixmap(_fromUtf8(saveas_img_path)),
                              QtGui.QIcon.Normal,
                              QtGui.QIcon.Off)
        self.saveFileAsAction = QtGui.QAction(icon_saveas, "Save as", self)
        self.saveFileAsAction.setStatusTip('Save file as')
        self.saveFileAsAction.triggered.connect(self.saveFileAs)


        # Close file action
        close_img_path = ":/icons/22x22/process-stop.png"
        icon_close     = QtGui.QIcon()
        icon_close.addPixmap(QtGui.QPixmap(_fromUtf8(close_img_path)),
                             QtGui.QIcon.Normal,
                             QtGui.QIcon.Off)
        self.closeFileAction  = QtGui.QAction(icon_close, "Close file", self)
        self.closeFileAction.setShortcut("Ctrl+Q")
        self.closeFileAction.setStatusTip('Close opened file')
        self.closeFileAction.triggered.connect(self.closeOpenedFile)

        # Exit editor action
        quit_img_path = ":/icons/22x22/system-log-out.png"
        icon_quit     = QtGui.QIcon()
        icon_quit.addPixmap(QtGui.QPixmap(_fromUtf8(quit_img_path)),
                          QtGui.QIcon.Normal,
                          QtGui.QIcon.Off)
        self.quitAction = QtGui.QAction(icon_quit, "Quit", self)
        self.quitAction.setStatusTip('Quit the editor')
        self.quitAction.triggered.connect(self.closeApplication)

        self.statusBar()

        # File toolbar
        self.toolbar = self.addToolBar("Options")

        self.toolbar.addAction(self.newFileAction)
        self.toolbar.addAction(self.openFileAction)
        self.toolbar.addAction(self.saveFileAction)
        self.toolbar.addAction(self.saveFileAsAction)
        self.toolbar.addAction(self.closeFileAction)
        self.toolbar.addAction(self.quitAction)

        # File menu
        self.mainMenu = self.menuBar()

        self.fileMenu = self.mainMenu.addMenu('&File')
        self.fileMenu.addAction(self.newFileAction)
        self.fileMenu.addAction(self.openFileAction)
        self.fileMenu.addAction(self.saveFileAction)
        self.fileMenu.addAction(self.saveFileAsAction)
        self.fileMenu.addAction(self.closeFileAction)
        self.fileMenu.addAction(self.quitAction)

        # Explorer
        self.explorer = self._initFileExplorer()
        self._initExplorerActions()

        # Editor
        self.textEdit = self._initFileEditor()

        # file attributes
        self.filename = ""

        self.mainWidget = FormWidget(self, [self.explorer, self.textEdit])
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
        textEdit = QtGui.QTextEdit()
        textEdit.setFont(base_font)
        textEdit.textChanged.connect(self.updateFileState)
        textEdit.setReadOnly(self.readOnly)
        policy = textEdit.sizePolicy()
        policy.setHorizontalPolicy(QtGui.QSizePolicy.Expanding)
        textEdit.setSizePolicy(policy)

        # tab
        textEdit.setTabStopWidth(font_metrics.width(_tab_string))

        return textEdit
    # ---------------------------------------------------------------


    # ---------------------------------------------------------------
    def _initFileExplorer(self, base_dir=None):
        """
        Create the File explorer object based on the QFileSystemModel widget.
        """

        model = QtWidgets.QFileSystemModel()
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

        return tree;
    # ---------------------------------------------------------------


    # ---------------------------------------------------------------
    def _initExplorerActions(self):
        """
        Create explorer actions dictionary
        """

        _editAction = QtGui.QAction(self.explorer.model())
        _editAction.setText('Edit file')
        _editAction.triggered.connect(self._editSelectedFile)

        _viewAction = QtGui.QAction(self.explorer.model())
        _viewAction.setText('View file')
        _viewAction.triggered.connect(self._viewSelectedFile)

        _copyAction = QtGui.QAction(self.explorer.model())
        _copyAction.setText('Copy to SRC')
        _copyAction.triggered.connect(self._copySelectedFile)

        _deleteAction = QtGui.QAction(self.explorer.model())
        _deleteAction.setText('Remove from SRC')
        _deleteAction.triggered.connect(self._removeSelectedFile)

        self._explorerActions = {'edit':_editAction,
                                 'view':_viewAction,
                                 'copy':_copyAction,
                                 'remove':_deleteAction}
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

        if self.case_name == 'SRC':
            trg_path = os.path.join(self.case_dir,
                                    self._currentSelection['filename'])
        else:
            sp = self._currentSelection['subpath']
            while '/' in sp and len(sp) > 3:
                e1, e2 = os.path.split(sp)
                if e2 == 'SRC':
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
        question = "Remove %s from the SRC folder ?" % (self._currentSelection['filename'])

        choice = QtGui.QMessageBox.question(self,
                                            title,
                                            question,
                                            QtGui.QMessageBox.Yes | QtGui.QMessageBox.No)

        if choice == QtGui.QMessageBox.Yes:
            fn = os.path.join(self.case_dir,
                              self._currentSelection['subpath'],
                              self._currentSelection['filename'])

            os.remove(fn)
        else:
            pass
    # ---------------------------------------------------------------


    # ---------------------------------------------------------------
    def explorerContextMenu(self, position):
        """
        Custom menu for the mouse right-click.
        Depends on whether the file is in the SRC, SRC/subfolder
        or RESU/subfolder.
        Possible actions are 'edit', 'view' and 'copy' (to SRC)
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
                if ptxt in [None, self.case_name]:
                    pe = ptxt
                    break
                path2file = ptxt + '/' + path2file
                c = p
                p = c.parent()

        self._currentSelection = {'filename':fname,
                                  'subpath' :path2file,
                                  'filedir' :ps,
                                  'origdir' :pe}


        self._contextMenu = QtGui.QMenu()
        if pe == 'RESU':
            self._contextMenu.addAction(self._explorerActions['view'])
        elif pe == 'SRC':
            if ps == 'SRC':
                self._contextMenu.addAction(self._explorerActions['edit'])
                self._contextMenu.addAction(self._explorerActions['remove'])
            else:
                self._contextMenu.addAction(self._explorerActions['view'])
                self._contextMenu.addAction(self._explorerActions['copy'])

        self._contextMenu.exec_(self.explorer.viewport().mapToGlobal(position))
    # ---------------------------------------------------------------


    # ---------------------------------------------------------------
    def updateFileState(self, new_state = False):
        """
        Update file state (saved or not)
        """
        self.saved = new_state
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
            self.filename = self.dialog.getOpenFileName(self, 'Open File', self.last_dir)

        self.last_dir = os.path.split(self.filename)[0]

        self.textEdit.setReadOnly(self.readOnly)
        self.saveFileAction.setEnabled(not self.readOnly)
        self.saveFileAsAction.setEnabled(not self.readOnly)

        if self.filename != None and self.filename != '':
            file = open(self.filename,'r')

            self.newFile()
            with file:
                text = file.read()
                self.textEdit.setText(text)
                self.updateFileState(True)
    # ---------------------------------------------------------------


    # ---------------------------------------------------------------
    def newFile(self):
        """
        Create a new file (blank)
        """

        self.updateFileState(False)
        hl = QtextHighlighter(self.textEdit)
        self.textEdit.show()
    # ---------------------------------------------------------------


    # ---------------------------------------------------------------
    def saveFile(self):
        """
        Save file
        """

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

        self.filename = self.dialog.getSaveFileName(self, 'Save File')
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

        if self.saved == False:
            choice = QtGui.QMessageBox.question(self, 'Built-in editor',
                                                'File changed.\nDo you want to save?',
                                                QtGui.QMessageBox.Yes | QtGui.QMessageBox.No)
            if choice == QtGui.QMessageBox.Yes:
                self.saveFile()
            else:
                pass

        self.saved = True
        self.filename = ''
        self.textEdit.setText('')
    # ---------------------------------------------------------------


    # ---------------------------------------------------------------
    def closeApplication(self):
        """
        Close the editor
        """
        choice = QtGui.QMessageBox.question(self, 'Built-in editor',
                                            "Exit text editor?",
                                            QtGui.QMessageBox.Yes | QtGui.QMessageBox.No)
        if choice == QtGui.QMessageBox.Yes:
            self.closeOpenedFile()
            self.close()
        else:
            pass
    # ---------------------------------------------------------------


#-------------------------------------------------------------------------------
# END OF FILE
#-------------------------------------------------------------------------------
