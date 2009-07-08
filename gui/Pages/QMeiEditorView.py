# -*- coding: utf-8 -*-
#
#-------------------------------------------------------------------------------
#
#     This file is part of the Code_Saturne User Interface, element of the
#     Code_Saturne CFD tool.
#
#     Copyright (C) 1998-2007 EDF S.A., France
#
#     contact: saturne-support@edf.fr
#
#     The Code_Saturne User Interface is free software; you can redistribute it
#     and/or modify it under the terms of the GNU General Public License
#     as published by the Free Software Foundation; either version 2 of
#     the License, or (at your option) any later version.
#
#     The Code_Saturne User Interface is distributed in the hope that it will be
#     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
#     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with the Code_Saturne Kernel; if not, write to the
#     Free Software Foundation, Inc.,
#     51 Franklin St, Fifth Floor,
#     Boston, MA  02110-1301  USA
#
#-------------------------------------------------------------------------------

"""
This module defines the Dialog window of the mathematical expression interpretor.

This module contains the following classes and function:
- QMeiEditorView
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import sys, string
import logging

#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------

from PyQt4.QtCore import *
from PyQt4.QtGui  import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from Base.Toolbox import GuiParam
from QMeiEditorForm import Ui_QMeiDialog
import mei

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("QMeiEditorView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# 
#-------------------------------------------------------------------------------

class QMeiHighlighter(QSyntaxHighlighter):

    def __init__(self, document, required, symbols):
        QSyntaxHighlighter.__init__(self, document)

        list = []
        for p,q in required:
            list.append(p)
        rx = string.join(list, r"|")
        self.req_re  = QRegExp(rx)
        self.req_fmt = QTextCharFormat()
        self.req_fmt.setFontWeight(QFont.Bold)

        list = []
        for p,q in symbols:
            list.append(p)
        rx = string.join(list, r"|")
        self.sym_re  = QRegExp(rx)
        self.sym_fmt = QTextCharFormat()
        self.sym_fmt.setFontWeight(QFont.Bold)

        self.const_re  = QRegExp(r"e|pi")
        self.const_fmt = QTextCharFormat()
        self.const_fmt.setForeground(Qt.darkCyan)

        self.func_re  = QRegExp(r"exp|log|sqrt|sin|cos|tan|asin|acos|atan|atan2|sinh|cosh|tanh|abs|min|max")
        self.func_fmt = QTextCharFormat()
        self.func_fmt.setFontWeight(QFont.Bold)
        self.func_fmt.setForeground(Qt.blue)

        self.kword_re  = QRegExp(r"while|if|else|print")
        self.kword_fmt = QTextCharFormat()
        self.kword_fmt.setForeground(Qt.darkRed)

        #self.special_re = QRegExp(r"[^\\][(),]")
        #self.special_fmt = QTextCharFormat()
        #self.special_fmt.setForeground(Qt.blue)

        self.rules = [
            (self.req_re,   self.req_fmt, 0, 0),
            (self.sym_re,   self.sym_fmt, 0, 0),
            (self.const_re, self.const_fmt, 0, 0),
            (self.func_re,  self.func_fmt, 0, 0),
            (self.kword_re, self.kword_fmt, 0, 0),
            #(self.special_re, self.special_fmt, 1, -1),
        ]


    def highlightBlock(self, text):
        for expr, fmt, a, b in self.rules:
            index = text.indexOf(expr)
            while index >= 0:
                length = expr.matchedLength()
                self.setFormat(index + a, length + b, fmt)
                index = text.indexOf(expr, index + length + b)

#-------------------------------------------------------------------------------
# 
#-------------------------------------------------------------------------------

class QMeiEditorView(QDialog, Ui_QMeiDialog):
    """
    """
    def __init__(self, parent, expression = "", symbols = [], required = [], examples = ""):
        """
        Constructor.
        """
        QDialog.__init__(self, parent)

        Ui_QMeiDialog.__init__(self)
        self.setupUi(self)

        self.required = required
        self.symbols  = symbols

        #h1 = QMeiHighlighter(self.textEditExpression, required, symbols)
        #h2 = QMeiHighlighter(self.textEditExamples, required, symbols)

        # Required symbols of the mathematical expression

        if not expression:
            expression = ""
            for p,q in required:
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
                  "<big><u>Usefull functions:</u></big><br>"\
                  "<b>cos</b>: cosine<br>"\
                  "<b>sin</b>: sinus<br>"\
                  "<b>tan</b>: tangent<br>"\
                  "<b>exp</b>: exponential<br>"\
                  "<b>sqrt</b>: square root<br>"\
                  "<b>log</b>: napierian logarithm<br>"\
                  "<b>acos</b>: arccosine<br>"\
                  "<b>asin</b>: arcsine<br>"\
                  "<b>atan</b>: arctangent<br>"\
                  "<b>atan2</b>: arctangent<br>"\
                  "<b>cosh</b>: hyperbolic cosine<br>"\
                  "<b>sinh</b>: hyperbolic sine<br>"\
                  "<b>tanh</b>: hyperbolic tangent<br>"\
                  "<b>abs</b>: absolute value<br>"\
                  "<b>mod</b>: modulo<br>"\
                  "<b>int</b>: floor<br>"\
                  "<b>min</b>: minimum<br>"\
                  "<b>max</b>: maximum<br>"\
                  "<br>"\
                  "<big><u>Usefull constants:</u></big><br>"\
                  "<b>pi</b> = 3.14159265358979323846<br>"\
                  "<b>e</b> = 2.718281828459045235<br>"\
                  "<br>"\
                  "<big><u>Operators and statements:</u></big><br>"\
                  "<b><code> + - * / ^ </code></b><br>"\
                  "<b><code>! &lt; &gt; &lt;= &gt;= == != && || </code></b><br>"\
                  "<b><code>while if else</code></b><br>"\
                  "")

        # lay out the text

        self.textEditExpression.setText(expression)
        self.textEditSymbols.setText(predif)
        self.textEditExamples.setText(examples)

        self.expressionDoc = self.textEditExpression.document()


    @pyqtSignature("")
    def slotClearBackground(self):
        """
        Private slot.
        """
        QObject.disconnect(self.textEditExpression, SIGNAL("textChanged()"), self.slotClearBackground)
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
        log.debug("accept()")

        doc = self.textEditExpression.document()

        intr = mei.mei_tree_new(str(self.textEditExpression.toPlainText()))
        log.debug("intr.string: %s" % intr.string)

        for p,q in self.symbols:
            mei.mei_tree_insert(intr, p, 0.0)

        # Unknown symbols

        if mei.mei_tree_builder(intr):
            log.debug("intr.errors: %s" % intr.errors)
            for i in range(0, intr.errors):
                log.debug("  line:   %s \n  column: %s\n" % (intr.lines[i], intr.columns[i]))
            msg = ""
            for i in range(0, intr.errors):
                l = intr.lines[i]
                c = intr.columns[i]
                block_text = doc.findBlockByNumber(l - 1)
                cursor = QTextCursor(block_text)
                block_format = QTextBlockFormat()
                block_format.setBackground(QBrush(QColor(Qt.red)))
                cursor.setBlockFormat(block_format)

                #self.textEditExpression.setCursorPosition(l, c)
                #self.textEditExpression.setFocus()
                msg += intr.labels[i] + \
                       "    line: "   + str(l)   + " \n" + \
                       "    column: " + str(c) + " \n\n"

            QMessageBox.critical(self, self.tr('Expression Editor'), QString(msg))
            intr = mei.mei_tree_destroy(intr)
            QObject.connect(self.textEditExpression, SIGNAL("textChanged()"), self.slotClearBackground)
            return

        # If required symbols are missing

        list = []
        for p,q in self.required:
            list.append(p)

        if mei.mei_tree_find_symbols(intr, len(list), list):
            msg = "Warning, these required symbols are not found:\n\n"
            for i in range(0, doc.blockCount()):
                block_text = doc.findBlockByNumber(i)
                cursor = QTextCursor(block_text)
                block_format = QTextBlockFormat()
                block_format.setBackground(QBrush(QColor(Qt.magenta)))
                cursor.setBlockFormat(block_format)
                #self.textEditExpression.setFocus()
            for i in range(0, intr.errors): msg += intr.labels[i] + " \n"
            QMessageBox.critical(self, self.tr('Expression Editor'), QString(msg))
            intr = mei.mei_tree_destroy(intr)
            QObject.connect(self.textEditExpression, SIGNAL("textChanged()"), self.slotClearBackground)
            return

        intr = mei.mei_tree_destroy(intr)

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
    app.connect(app, SIGNAL("lastWindowClosed()"), app, SLOT("quit()"))
    parent = QWidget()
    dlg = QMeiEditorView(parent)
    dlg.show()
    signal.signal(signal.SIGINT, signal.SIG_DFL)
    sys.exit(app.exec_())

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
