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

#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------

from code_saturne.gui.base.QtCore    import *
from code_saturne.gui.base.QtGui     import *
from code_saturne.gui.base.QtWidgets import *

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
        _format.setFontWeight(QFont.Bold)
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
        '!', '==', '!=', '<', '<=', '>', '>=', '&&', '\|\|',
        # Arithmetic
        '=', '\+', '-', '\*', '/', '\^',
    ]

    braces = ['\{', '\}', '\(', '\)', '\[', '\]',]


    def __init__(self, document, symbols):
        QSyntaxHighlighter.__init__(self, document)

        keywordFormat = QTextCharFormat()
        keywordFormat.setForeground(Qt.darkMagenta)
        keywordFormat.setFontWeight(QFont.Bold)

        keywordPatterns = ["\\b?xml\\b", "/>", ">", "<"]

        self.highlightingRules = [(QRegExp(pattern), keywordFormat)
                for pattern in keywordPatterns]

        xmlElementFormat = QTextCharFormat()
        xmlElementFormat.setFontWeight(QFont.Bold)
        xmlElementFormat.setForeground(Qt.green)
        self.highlightingRules.append((QRegExp("\\b[A-Za-z0-9_]+(?=[\s/>])"), xmlElementFormat))

        xmlAttributeFormat = QTextCharFormat()
        xmlAttributeFormat.setFontItalic(True)
        xmlAttributeFormat.setForeground(Qt.blue)
        self.highlightingRules.append((QRegExp("\\b[A-Za-z0-9_]+(?=\\=)"), xmlAttributeFormat))

        self.valueFormat = QTextCharFormat()
        self.valueFormat.setForeground(Qt.red)

        self.valueStartExpression = QRegExp("\"")
        self.valueEndExpression = QRegExp("\"(?=[\s></])")

        singleLineCommentFormat = QTextCharFormat()
        singleLineCommentFormat.setForeground(Qt.gray)
        self.highlightingRules.append((QRegExp("<!--[^\n]*-->"), singleLineCommentFormat))


    def highlightBlock(self, text):
        """
        Apply syntax highlighting to the given block of text.
        """
        for pattern, format in self.highlightingRules:

            #Create a regular expression from the retrieved pattern
            expression = QRegExp(pattern)

            #Check what index that expression occurs at with the ENTIRE text
            index = expression.indexIn(text)

            while index >= 0:
                length = expression.matchedLength()
                self.setFormat(index, length, format)

                index = expression.indexIn(text, index + length)

        self.setCurrentBlockState(0)

        startIndex = 0
        if self.previousBlockState() != 1:
            startIndex = self.valueStartExpression.indexIn(text)

        while startIndex >= 0:
            endIndex = self.valueEndExpression.indexIn(text, startIndex)

            if endIndex == -1:
                self.setCurrentBlockState(1)
                commentLength = len(text) - startIndex
            else:
                commentLength = endIndex - startIndex + self.valueEndExpression.matchedLength()

            self.setFormat(startIndex, commentLength, self.valueFormat)

            startIndex = self.valueStartExpression.indexIn(text, startIndex + commentLength);


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


    def accept(self):
        """
        What to do when user clicks on 'OK'.
        """
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




