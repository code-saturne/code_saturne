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
This module defines the Dialog window of the mathematical expression interpretor.

This module contains the following classes and function:
- format
- QMeiHighlighter
- QMeiEditorView
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import os
import sys
import logging
import subprocess

#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------

from code_saturne.Base.QtCore    import *
from code_saturne.Base.QtGui     import *
from code_saturne.Base.QtWidgets import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.Base.Toolbox import GuiParam
from code_saturne.Pages.QMeiEditorForm import Ui_QMeiDialog

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("QMeiEditorView")
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
    'required': format('magenta', 'bold'),
    'symbols': format('darkMagenta', 'bold'),
    'comment': format('darkGreen', 'italic'),
    }


class QMeiHighlighter(QSyntaxHighlighter):
    """
    Syntax highlighter for the mathematical expressions editor.
    """
    keywords = [
         'cos', 'sin', 'tan', 'exp', 'sqrt', 'log',
         'acos', 'asin', 'atan', 'atan2', 'cosh', 'sinh',
         'tanh', 'abs', 'mod', 'int', 'min', 'max',
         'pi', 'e', 'while', 'if', 'else', 'print',
         'while'
    ]

    operators = [
        # logical
        '!', '==', '!=', '<', '<=', '>', '>=', '&&', '\|\|',
        # Arithmetic
        '=', '\+', '-', '\*', '/', '\^',
    ]

    braces = ['\{', '\}', '\(', '\)', '\[', '\]',]


    def __init__(self, document, required, symbols):
        QSyntaxHighlighter.__init__(self, document)

        rules = []
        rules += [(r'\b%s\b' % w, STYLES['keyword']) for w in QMeiHighlighter.keywords]
        rules += [(r'%s' % o, STYLES['operator']) for o in QMeiHighlighter.operators]
        rules += [(r'%s' % b, STYLES['brace']) for b in QMeiHighlighter.braces]
        rules += [(r'\b%s\b' % r, STYLES['required']) for r,l in required]
        rules += [(r'\b%s\b' % s, STYLES['symbols']) for s,l in symbols]
        rules += [(r'#[^\n]*', STYLES['comment'])]

        # Build a QRegExp for each pattern
        self.rules = [(QRegExp(pat), fmt) for (pat, fmt) in rules]


    def highlightBlock(self, text):
        """
        Apply syntax highlighting to the given block of text.
        """
        for rx, fmt in self.rules:
            pos = rx.indexIn(text, 0)
            while pos != -1:
                pos = rx.pos(0)
                s = rx.cap(0)
                self.setFormat(pos, len(s), fmt)
                pos = rx.indexIn( text, pos+rx.matchedLength() )

#-------------------------------------------------------------------------------
# Dialog for mathematical expression interpretor
#-------------------------------------------------------------------------------

class QMeiEditorView(QDialog, Ui_QMeiDialog):
    """
    """
    def __init__(self, parent, check_syntax = None, expression = "", symbols = [], required = [], examples = ""):
        """
        Constructor.
        """
        QDialog.__init__(self, parent)

        Ui_QMeiDialog.__init__(self)
        self.setupUi(self)

        if check_syntax is not None:
            if not os.path.isfile(check_syntax):
                check_syntax = None

        self.check_syntax = check_syntax

        sym = [('todelete','')]
        base = []
        for s,l in symbols:
            if s.find("[") != -1:
                idx = s.find("[")
                if s[:idx] not in base:
                    sym.append((s[:idx], l))
                    base.append(s[:idx])
            else:
                sym.append((s,l))
        del sym[0]

        req = [('todelete','')]
        base = []
        for s,l in required:
            if s.find("[") != -1:
                idx = s.find("[")
                if s[:idx] not in base:
                    req.append((s[:idx], l))
                    base.append(s[:idx])
            else:
                req.append((s,l))
        del req[0]

        self.required = required
        self.symbols  = symbols

        # Syntax highlighting
        self.h1 = QMeiHighlighter(self.textEditExpression, req, sym)
        self.h2 = QMeiHighlighter(self.textEditExamples, req, sym)

        # Required symbols of the mathematical expression

        if not expression:
            expression = ""
            for p, q in required:
                expression += p + " = \n"
        else:
            while expression[0] in (" ", "\n", "\t"):
                expression = expression[1:]
            while expression[-1] in (" ", "\n", "\t"):
                expression = expression[:-1]

        # Predefined symbols for the mathematical expression

        predif = self.tr("<big><u>Required symbol:</u></big><br>")

        for p,q in required:
            for s in ("<b>", p, "</b>: ", q, "<br>"):
                predif += s

        if symbols:
            predif += self.tr("<br><big><u>Predefined symbols:</u></big><br>")

            for p,q in symbols:
                for s in ("<b>", p, "</b>: ", q, "<br>"):
                    predif += s

        predif += self.tr("<br>"\
                  "<big><u>Useful functions:</u></big><br>"\
                  "<b>cos</b>: cosine<br>"\
                  "<b>sin</b>: sine<br>"\
                  "<b>tan</b>: tangent<br>"\
                  "<b>exp</b>: exponential<br>"\
                  "<b>sqrt</b>: square root<br>"\
                  "<b>log</b>: napierian logarithm<br>"\
                  "<b>acos</b>: arc cosine<br>"\
                  "<b>asin</b>: arcsine<br>"\
                  "<b>atan</b>: arc tangent<br>"\
                  "<b>atan2</b>: arc tangent (two variables)<br>"\
                  "<b>cosh</b>: hyperbolic cosine<br>"\
                  "<b>sinh</b>: hyperbolic sine<br>"\
                  "<b>tanh</b>: hyperbolic tangent<br>"\
                  "<b>abs</b>: absolute value<br>"\
                  "<b>mod</b>: modulo<br>"\
                  "<b>int</b>: floor<br>"\
                  "<b>min</b>: minimum<br>"\
                  "<b>max</b>: maximum<br>"\
                  "<br>"\
                  "<big><u>Useful constants:</u></big><br>"\
                  "<b>pi</b> = 3.14159265358979323846<br>"\
                  "<b>e</b> = 2.718281828459045235<br>"\
                  "<br>"\
                  "<big><u>Operators and statements:</u></big><br>"\
                  "<b><code> + - * / ^ </code></b><br>"\
                  "<b><code>! &lt; &gt; &lt;= &gt;= == != && || </code></b><br>"\
                  "<b><code>while if else print</code></b><br>"\
                  "")

        # lay out the text

        self.textEditExpression.setText(expression)
        self.textEditSymbols.setText(predif)
        self.textEditExamples.setText(examples)

        self.expressionDoc = self.textEditExpression.document()


    @pyqtSlot()
    def slotClearBackground(self):
        """
        Private slot.
        """
        QObject.disconnect(self.textEditExpression)
        doc = self.textEditExpression.document()

        for i in range(0, doc.blockCount()):
            block_text = doc.findBlockByNumber(i)
            log.debug("block: %s" % block_text.text())
            cursor = QTextCursor(block_text)
            block_format = QTextBlockFormat()
            block_format.clearBackground()
            cursor.setBlockFormat(block_format)


    def accept(self):
        """
        What to do when user clicks on 'OK'.
        """

        if self.check_syntax == None:
            QDialog.accept(self)
            return

        log.debug("accept()")

        doc = self.textEditExpression.document()

        log.debug("check.string: %s" % str(self.textEditExpression.toPlainText()))

        check = subprocess.Popen([self.check_syntax],
                                 stdin=subprocess.PIPE,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 universal_newlines=True)
        check.stdin.write(str(self.textEditExpression.toPlainText()) + '\30')

        check.stdin.write(str(len(self.symbols)) + '\30')
        for p,q in self.symbols:
            check.stdin.write(str(p) + '\30')

        check.stdin.write(str(len(self.required)) + '\30')
        for p,q in self.required:
            check.stdin.write(str(p) + '\30')

        check_out, check_err = check.communicate()

        # Invalid expression

        if check.returncode == 2:

            # In this case, first record defines the number of errors,
            # and following records given by groups of 3,
            # for line, column, and label.

            errors = check_err.split('\30')
            log.debug("check.errors: %s" % errors)
            n_errors = int(errors[0])

            msg = ""
            for i in range(0, n_errors):
                l = int(errors[i*3+1])
                c = int(errors[i*3+2])
                log.debug("  line:   %s \n  column: %s\n" % (l, c))
                block_text = doc.findBlockByNumber(l - 1)
                cursor = QTextCursor(block_text)
                block_format = QTextBlockFormat()
                block_format.setBackground(QBrush(QColor(Qt.red)))
                cursor.setBlockFormat(block_format)

                msg += errors[i*3+3] + \
                       "    line: "   + str(l)   + " \n" + \
                       "    column: " + str(c) + " \n\n"

            QMessageBox.critical(self, self.tr('Expression Editor'), str(msg))
            self.textEditExpression.textChanged.connect(self.slotClearBackground)
            return

        # If required symbols are missing

        elif check.returncode == 3:

            # In this case, first record defines the number of errors,
            # and following records give matching labels.

            errors = check_err.split('\30')
            log.debug("check.errors: %s" % errors)
            n_errors = int(errors[0])

            msg = "Warning, these required symbols are not found:\n\n"
            for i in range(0, doc.blockCount()):
                block_text = doc.findBlockByNumber(i)
                cursor = QTextCursor(block_text)
                block_format = QTextBlockFormat()
                block_format.setBackground(QBrush(QColor(Qt.yellow)))
                cursor.setBlockFormat(block_format)
            for i in range(0, n_errors): msg += errors[i+1] + " \n"
            QMessageBox.critical(self, self.tr('Expression Editor'), str(msg))
            self.textEditExpression.textChanged.connect(self.slotClearBackground)
            return

        # If another error was found

        elif check.returncode != 0:

            errors = check_err.split('\30')
            log.debug("check.errors: %s" % errors)

            msg = "Warning, expression check failed unexpectedly:\n\n"
            QMessageBox.critical(self, self.tr('Expression Editor'), str(msg))
            self.textEditExpression.textChanged.connect(self.slotClearBackground)
            return

        QDialog.accept(self)


    def reject(self):
        """
        Method called when 'Cancel' button is clicked
        """
        log.debug("reject()")
        QDialog.reject(self)


    def get_result(self):
        """
        Method to get the result
        """
        return self.textEditExpression.toPlainText()


    def tr(self, text):
        """
        Translation
        """
        return text

#-------------------------------------------------------------------------------
# Test function
#-------------------------------------------------------------------------------

if __name__ == "__main__":
    import sys, signal
    app = QApplication(sys.argv)
    app.lastWindowClosed.connect(app.quit)
    parent = QWidget()
    dlg = QMeiEditorView(parent)
    dlg.show()
    signal.signal(signal.SIGINT, signal.SIG_DFL)
    sys.exit(app.exec_())

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
