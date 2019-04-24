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

import sys
from code_saturne.Base import QtGui, QtCore

# Check if QString exists
has_qstring = True
try:
    from code_saturne.Base.QtCore import QString
    _fromUtf8 = QString.fromUtf8
except ImportError:
    has_qstring = False
    def _fromUtf8(s):
        return s

import resource_base_rc
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

class HighlightingRule():

    def __init__(self, pattern, format):

        self.pattern = pattern
        self.format  = format
#-------------------------------------------------------------------------------
# QtextHighlighter class
#-------------------------------------------------------------------------------

class QtextHighlighter(QtGui.QSyntaxHighlighter):
    """
    Syntax highighting
    """

    def __init__(self, parent):

        QtGui.QSyntaxHighlighter.__init__(self, parent)
        self.parent = parent
        self.highlightingRules = []

        # Keywords (C or Fortran)
        self.kw = ['if', 'else', 'endif', '\#', 'include',
                   'void', 'int', 'integer', 'double', 'const',
                   'fprintf', 'bft_printf', 'bft_printf_flush',
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


#-------------------------------------------------------------------------------
# QFileEditor class
#-------------------------------------------------------------------------------

class QFileEditor(QtGui.QMainWindow):

    def __init__(self, parent=None):
        super(QFileEditor, self).__init__(parent)
        self.setGeometry(50, 50, 500, 300)
        self.setWindowTitle("Text editor")

        # Open file action
        open_img_path = ":/icons/22x22/document-open.png"
        icon_open     = QtGui.QIcon()
        icon_open.addPixmap(QtGui.QPixmap(_fromUtf8(open_img_path)),
                            QtGui.QIcon.Normal,
                            QtGui.QIcon.Off)
        openFile = QtGui.QAction(icon_open, "Open", self)
        openFile.setShortcut("Ctrl+O")
        openFile.setStatusTip('Open File')
        openFile.triggered.connect(self.openFile)

        # New file action
        new_img_path = ":/icons/22x22/document-new.png"
        icon_new     = QtGui.QIcon()
        icon_new.addPixmap(QtGui.QPixmap(_fromUtf8(new_img_path)),
                          QtGui.QIcon.Normal,
                          QtGui.QIcon.Off)
        newFile = QtGui.QAction(icon_new, "New", self)
        newFile.setShortcut("Ctrl+E")
        newFile.setStatusTip('Create new file')
        newFile.triggered.connect(self.newFile)

        # Save as action
        save_img_path = ":/icons/22x22/document-save-as.png"
        icon_save     = QtGui.QIcon()
        icon_save.addPixmap(QtGui.QPixmap(_fromUtf8(save_img_path)),
                          QtGui.QIcon.Normal,
                          QtGui.QIcon.Off)
        saveFile = QtGui.QAction(icon_save, "Save", self)
        saveFile.setShortcut("Ctrl+S")
        saveFile.setStatusTip('Save file')
        saveFile.triggered.connect(self.saveFile)

        # Exit editor action
        quit_img_path = ":/icons/22x22/system-log-out.png"
        icon_quit     = QtGui.QIcon()
        icon_quit.addPixmap(QtGui.QPixmap(_fromUtf8(quit_img_path)),
                          QtGui.QIcon.Normal,
                          QtGui.QIcon.Off)
        quitAction = QtGui.QAction(icon_quit, "Quit", self)
        quitAction.setShortcut("Ctrl+Q")
        quitAction.setStatusTip('Quit the editor')
        quitAction.triggered.connect(self.closeApplication)

        self.statusBar()

        # File toolbar
        self.toolbar = self.addToolBar("Options")

        self.toolbar.addAction(newFile)
        self.toolbar.addAction(openFile)
        self.toolbar.addAction(saveFile)
        self.toolbar.addAction(quitAction)

        # File menu
        mainMenu = self.menuBar()

        fileMenu = mainMenu.addMenu('&File')
        fileMenu.addAction(newFile)
        fileMenu.addAction(openFile)
        fileMenu.addAction(saveFile)
        fileMenu.addAction(quitAction)

        self.textEdit = QtGui.QTextEdit()
        self.textEdit.textChanged.connect(self.updateFileState)

    def updateFileState(self, new_state = False):
        self.saved = new_state

    def openFile(self):
        name = QtGui.QFileDialog.getOpenFileName(self, 'Open File')

        if name != None and name != '':
            file = open(name,'r')

            self.newFile()
            with file:
                text = file.read()
                self.textEdit.setText(text)
                self.saved = True


    def newFile(self):

        self.saved = False
        hl = QtextHighlighter(self.textEdit)
        self.setCentralWidget(self.textEdit)
        self.textEdit.show()



    def saveFile(self):
        name = QtGui.QFileDialog.getSaveFileName(self, 'Save File')

        if name != None and name != '':
            file = open(name,'w')
            text = self.textEdit.toPlainText()
            file.write(text)
            file.close()

            self.saved = True


    def closeApplication(self):
        choice = QtGui.QMessageBox.question(self, 'Built-in editor',
                                            "Exit text editor?",
                                            QtGui.QMessageBox.Yes | QtGui.QMessageBox.No)
        if choice == QtGui.QMessageBox.Yes:
#            sys.exit()
            if self.saved == False:
                choice = QtGui.QMessageBox.question(self, 'Built-in editor',
                                                    'File changed.\nDo you want to save?',
                                                    QtGui.QMessageBox.Yes | QtGui.QMessageBox.No)
                if choice == QtGui.QMessageBox.Yes:
                    self.saveFile()
                else:
                    pass
            self.close()
        else:
            pass


