"""
Provides QtCore classes and functions.
"""
import os
QT_API = "PYQT4"
try:
    from PyQt5.QtCore import QT_VERSION_STR
    QT_API = "PYQT5"
except:
    pass

if QT_API == "PYQT5":
    from PyQt5.QtCore import *
elif QT_API == "PYQT4":
    from PyQt4.QtCore import *
