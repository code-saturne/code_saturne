"""
Provide a supercharged QTextEdit method to allow for Autocompletion
"""

from types import MethodType

from QtGui import *
from QtWidgets import *
from QtCore import *
# ------------------------------------------------------------------------------
# QTextEdit with autocompletion
def CompletionTextEdit(target):

    target.completer = None

    def setCompleter(target, completer):
        if target.completer:
            target.completer.activated.disconnect()
        if not completer:
            return

        completer.setWidget(target)
        completer.setCompletionMode(QCompleter.PopupCompletion)
        completer.setCaseSensitivity(Qt.CaseInsensitive)
        target.completer = completer
        completer.setWidget(target)
        completer.activated.connect(target.insertCompletion)

    def insertCompletion(target, completion):
        tc = target.textCursor()
        if QT_API == "PYQT4":
            extra = (completion.length() -
                target.completer.completionPrefix().length())
            tc.movePosition(QTextCursor.Left)
            tc.movePosition(QTextCursor.EndOfWord)
            tc.insertText(completion.right(extra))
            target.setTextCursor(tc)
        elif QT_API == "PYQT5":
            extra = (len(completion) -
                len(target.completer.completionPrefix()))
            tc.movePosition(QTextCursor.Left)
            tc.movePosition(QTextCursor.EndOfWord)
            tc.insertText(completion[-extra:])
            target.setTextCursor(tc)

    def textUnderCursor(target):
        tc = target.textCursor()
        tc.select(QTextCursor.WordUnderCursor)
        return tc.selectedText()

    def focusInEvent(target, event):
        if target.completer:
            target.completer.setWidget(target);
        QTextEdit.focusInEvent(target, event)

    def keyPressEvent(target, event):
        if target.completer and target.completer.popup().isVisible():
            if event.key() in (
            Qt.Key_Enter,
            Qt.Key_Return,
            Qt.Key_Escape,
            Qt.Key_Tab,
            Qt.Key_Backtab):
                event.ignore()
                return

        ## has ctrl-E been pressed??
        isShortcut = (event.modifiers() == Qt.ControlModifier and
                      event.key() == Qt.Key_E)
        if (not target.completer or not isShortcut):
            QTextEdit.keyPressEvent(target, event)

        ## ctrl or shift key on it's own??
        ctrlOrShift = event.modifiers() in (Qt.ControlModifier ,
                    Qt.ShiftModifier)
        if QT_API == "PYQT4":
            if ctrlOrShift and event.text().isEmpty():
                # ctrl or shift key on it's own
                return
        elif QT_API == "PYQT5":
            if ctrlOrShift and len(event.text()) < 1:
                # ctrl or shift key on it's own
                return


        hasModifier = ((event.modifiers() != Qt.NoModifier) and
                        not ctrlOrShift)

        completionPrefix = target.textUnderCursor()

        # EOW test and compatibily with PyQt4/PyQt5
        if QT_API == "PYQT4":
            eow = QString("~!@#$%^&*()_+{}|:\"<>?,./;'[]\\-=") #end of word
            if (not isShortcut and (hasModifier or event.text().isEmpty() or
            completionPrefix.length() < 2 or
            eow.contains(event.text().right(1)))):
                target.completer.popup().hide()
                return
        elif QT_API == "PYQT5":
            eow = "~!@#$%^&*()_+{}|:\"<>?,./;'[]\\-=" #end of word
            if (not isShortcut and (hasModifier or len(event.text()) < 1 or
            len(completionPrefix) < 2 or
            event.text()[-1] in eow)):
                target.completer.popup().hide()
                return

        if (completionPrefix != target.completer.completionPrefix()):
            target.completer.setCompletionPrefix(completionPrefix)
            popup = target.completer.popup()
            popup.setCurrentIndex(
                target.completer.completionModel().index(0,0))

        cr = target.cursorRect()
        cr.setWidth(target.completer.popup().sizeHintForColumn(0)
            + target.completer.popup().verticalScrollBar().sizeHint().width())
        target.completer.complete(cr) ## popup it up!

    # Define methods
    target.setCompleter = MethodType(setCompleter,target)
    target.insertCompletion = MethodType(insertCompletion,target)
    target.textUnderCursor = MethodType(textUnderCursor,target)
    target.focusInEvent = MethodType(focusInEvent,target)
    target.keyPressEvent = MethodType(keyPressEvent,target)
