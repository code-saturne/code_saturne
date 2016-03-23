"""
Provides QtWidgets classes and functions.
"""
import os
QT_API = "PYQT4"
try:
    from PyQt5.QtCore import QT_VERSION_STR
    QT_API = "PYQT5"
except:
    pass

if QT_API == "PYQT5":
    from PyQt5.QtWidgets import *
elif QT_API == "PYQT4":
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

