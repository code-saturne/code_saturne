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
- SearchBar
Adapted from: https://www.binpress.com/building-text-editor-pyqt-3/ (08 July 2022)
"""

import sys, os, shutil
import re
from code_saturne.gui.base.QtWidgets import *
from code_saturne.gui.base.QtCore    import *
from code_saturne.gui.base.QtGui     import *

class SearchBar(QWidget):

    def __init__(self, textEdit, parent=None):
        """
        Initialize search bar. The search bar works on the content of textEdit.toPlainText().
        """
        QWidget.__init__(self, parent)
        self.textEditContent = textEdit
        self.lastMatch = None
        self.initUI()

    def initUI(self):
        """
        Create widgets and connect them to slots
        """
        self.lineEditFind = QLineEdit(self)
        self.lineEditFind.returnPressed.connect(self.find)
        self.pushButtonFind = QPushButton("Find", self)
        self.checkBoxCaseSensitive = QCheckBox("case sensitive", self)
        self.checkBoxWholeWords = QCheckBox("whole words", self)
        self.labelOccurences = QLabel("", self)

        layout = QGridLayout(self)
        layout.addWidget(self.lineEditFind, 0, 0, 1, 1)
        layout.addWidget(self.pushButtonFind, 0, 1, 1, 1)
        layout.addWidget(self.checkBoxCaseSensitive, 0, 2, 1, 1)
        layout.addWidget(self.checkBoxWholeWords, 0, 3, 1, 1)
        layout.addWidget(self.labelOccurences, 1, 0, 1, -1)

        self.lastMatch = None
        self.pushButtonFind.clicked.connect(self.find)
        shortcut = QShortcut(QKeySequence("Ctrl+F"), self)
        shortcut.activated.connect(self.startSearch)

    def find(self):
        """
        Find pattern with re library
        """
        text = self.textEditContent.toPlainText()
        query = self.lineEditFind.text()

        # Search options
        if self.checkBoxWholeWords.isChecked():
            query = r'\W' + query + r'\W'

        flags = re.I
        if self.checkBoxCaseSensitive.isChecked():
            flags = 0

        # Compile the regex pattern
        pattern = re.compile(query, flags)

        start = self.lastMatch.start() + 1 if self.lastMatch else 0
        self.lastMatch = pattern.search(text, start)
        nb_occurences = len(pattern.findall(text))
        if nb_occurences > 1:
            self.labelOccurences.setText("{0} occurences found".format(nb_occurences))
        else:
            self.labelOccurences.setText("{0} occurence found".format(nb_occurences))

        if self.lastMatch:
            start = self.lastMatch.start()
            end = self.lastMatch.end()
            # Remove '\W' from regex pattern for display
            if self.checkBoxWholeWords.isChecked():
                start += 1
                end -= 1
            self.moveCursor(start, end) # Beware ! Custom move cursor method

    def moveCursor(self, start, end):
        cursor = self.textEditContent.textCursor()
        cursor.setPosition(start)

        cursor.movePosition(QTextCursor.Right, QTextCursor.KeepAnchor, end - start)
        self.textEditContent.setTextCursor(cursor)

    def startSearch(self):
        if self.lineEditFind.text() != "":
            self.lineEditFind.selectAll()
        else:
            self.lineEditFind.setFocus()

    def hasSearchFocus(self):
        return self.lineEditFind.hasFocus()


if __name__ == "__main__":
    app = QApplication(sys.argv)
    SearchBar = SearchBar(None)
    SearchBar.show()
    sys.exit(app.exec_())
