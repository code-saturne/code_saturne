"""
Provides QtWidgets classes and functions.
"""

QT_API = "PYQT4"

from PyQt4.QtGui import *
from PyQt4.QtGui import QFileDialog as QFileDialogQt4

class QFileDialog(QFileDialogQt4):

    @staticmethod
    def getOpenFileName(parent=None, caption='', directory='',
                        filter='', selectedFilter='',
                        options=QFileDialogQt4.Options()):
        return QFileDialogQt4.getOpenFileNameAndFilter(
            parent, caption, directory, filter, selectedFilter,
            options)

    @staticmethod
    def getOpenFileNames(parent=None, caption='', directory='',
                         filter='', selectedFilter='',
                         options=QFileDialogQt4.Options()):
        return QFileDialogQt4.getOpenFileNamesAndFilter(
            parent, caption, directory, filter, selectedFilter,
            options)

    @staticmethod
    def getSaveFileName(parent=None, caption='', directory='',
                        filter='', selectedFilter='',
                        options=QFileDialogQt4.Options()):
        return QFileDialogQt4.getSaveFileNameAndFilter(
            parent, caption, directory, filter, selectedFilter,
            options)
