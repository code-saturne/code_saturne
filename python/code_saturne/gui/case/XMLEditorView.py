# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2024 EDF S.A.
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
This module defines the Dialog window of the XML viewer

This module contains the following classes and function:
- XMLHighlighter
- XMLEditorView
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import sys
import logging
import re

#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------

from code_saturne.gui.base.QtCore    import *
from code_saturne.gui.base.QtGui     import *
from code_saturne.gui.base.QtWidgets import *
from code_saturne.gui.base.SearchBar import SearchBar

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.model.Common import GuiParam
from code_saturne.gui.case.XMLEditorForm import Ui_XMLEditor

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("XMLEditorView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Syntax highlighter for the mathematical expressions editor.
#-------------------------------------------------------------------------------

def format(color, style=''):
    """Return a QTextCharFormat with the given attributes."""
    _color = QColor()
    _color.setNamedColor(color)

    _format = QTextCharFormat()
    _format.setForeground(_color)
    if 'bold' in style:
        _format.setFontWeight(QFont.Weight.Bold)
    if 'italic' in style:
        _format.setFontItalic(True)

    return _format


   # Syntax styles that can be shared by all expressions
STYLES = {
    'keyword': format('blue', 'bold'),
    'operator': format('red'),
    'brace': format('darkGray'),
    'symbols': format('darkMagenta', 'bold'),
    'comment': format('darkGreen', 'italic'),
    'numbers': format('brown'),
    }


class XMLHighlighter(QSyntaxHighlighter):
    """
    Syntax highlighter for the mathematical expressions editor.
    """
    keywords = [
         'solution_domain', 'velocity_pressure', 'variable', 'property',
         'turbulence', 'numerical_parameters'
    ]

    operators = [
        # logical
        '!', '==', '!=', '<', '<=', '>', '>=', '&&', r'\|\|',
        # Arithmetic
        '=', r'\+', '-', r'\*', '/', r'\^',
    ]

    braces = [r'\{', r'\}', r'\(', r'\)', r'\[', r'\]',]


    def __init__(self, document, symbols):
        QSyntaxHighlighter.__init__(self, document)

        keywordFormat = QTextCharFormat()
        keywordFormat.setForeground(Qt.GlobalColor.darkMagenta)
        keywordFormat.setFontWeight(QFont.Weight.Bold)

        keywordPatterns = [r"\b\?xml\b", "/>", ">", "<"]

        self.highlightingRules = [(QRegularExpression(pattern), keywordFormat)
                for pattern in keywordPatterns]

        xmlElementFormat = QTextCharFormat()
        xmlElementFormat.setFontWeight(QFont.Weight.Bold)
        xmlElementFormat.setForeground(Qt.GlobalColor.green)
        self.highlightingRules.append((QRegularExpression(r"\b[A-Za-z0-9_]+(?=[\s/>])"), xmlElementFormat))

        xmlAttributeFormat = QTextCharFormat()
        xmlAttributeFormat.setFontItalic(True)
        xmlAttributeFormat.setForeground(Qt.GlobalColor.blue)
        self.highlightingRules.append((QRegularExpression(r"\b[A-Za-z0-9_]+(?=\=)"), xmlAttributeFormat))

        self.valueFormat = QTextCharFormat()
        self.valueFormat.setForeground(Qt.GlobalColor.red)

        self.valueStartExpression = QRegularExpression("\"")
        self.valueEndExpression = QRegularExpression("\"(?=[\\s></])")

        singleLineCommentFormat = QTextCharFormat()
        singleLineCommentFormat.setForeground(Qt.GlobalColor.gray)
        self.highlightingRules.append((QRegularExpression("<!--[^\n]*-->"), singleLineCommentFormat))


    def highlightBlock(self, text):
        """
        Apply syntax highlighting to the given block of text.
        """
        for pattern, format in self.highlightingRules:

            #Create a regular expression from the retrieved pattern
            expression = QRegularExpression(pattern)

            #Check what index that expression occurs at with the ENTIRE text
            i = expression.globalMatch(text)
            while i.hasNext():
                match = i.next()
                index = match.capturedStart()
                length = match.capturedLength()
                self.setFormat(index, length, format)

        self.setCurrentBlockState(0)

        startIndex = 0
        if self.previousBlockState() != 1:
            match = self.valueStartExpression.match(text)
            if match.hasMatch():
                startIndex = match.capturedStart(0)
            else:
                startIndex = -1

        while startIndex >= 0:
            match = self.valueEndExpression.match(text, startIndex)
            if match.hasMatch():
                endIndex = match.capturedStart()
                startIndexNext = endIndex
                commentLength = endIndex - startIndex + match.capturedLength(0)
            else:
                endIndex = len(text)
                commentLength = endIndex - startIndex

            self.setFormat(startIndex, commentLength, self.valueFormat)

            startIndex += commentLength
            match = self.valueStartExpression.match(text, startIndex)
            if match.hasMatch():
                startIndex = match.capturedStart(0)
            else:
                startIndex = -1

#-------------------------------------------------------------------------------
# Dialog to show current XML status
#-------------------------------------------------------------------------------

class XMLEditorView(QDialog, Ui_XMLEditor):
    """
    """
    def __init__(self, parent, case):
        """
        Constructor.
        """
        QDialog.__init__(self, parent)

        Ui_XMLEditor.__init__(self)
        self.setupUi(self)
        self.create_widgets()

        self.symbols  = []
        self.case = case

        title = self.tr("XML source")
        self.setWindowTitle(title)

        # Syntax highlighting
        self.h1 = XMLHighlighter(self.textEditContent, self.symbols)

        # Current XML

        expression = self.case.toPrettyString()

        self.pushButtonValidate.clicked.connect(self.accept)

        # lay out the text

        self.textEditContent.setText(expression)

        self.expressionDoc = self.textEditContent.document()

    def create_widgets(self):
        """
        Add widgets programmatically
        """
        self.searchBar = SearchBar(self.textEditContent)
        self.layout().addWidget(self.searchBar, 0, 0, 1, -1)

    def accept(self):

        if self.searchBar.hasSearchFocus():
            self.searchBar.find()
        else:
            QDialog.accept(self)
        return


#-------------------------------------------------------------------------------
# Test function
#-------------------------------------------------------------------------------

if __name__ == "__main__":
    import sys, signal
    app = QApplication(sys.argv)
    app.lastWindowClosed.connect(app.quit)
    parent = QWidget()
    dlg = XMLEditorView(parent)
    dlg.show()
    signal.signal(signal.SIGINT, signal.SIG_DFL)
    sys.exit(app.exec_())

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------




