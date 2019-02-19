"""
Provide a supercharged QTextEdit method to allow for Autocompletion
"""

from types import MethodType

from QtGui import *
from QtWidgets import *
import QtCore
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
        completer.setCaseSensitivity(QtCore.Qt.CaseInsensitive)
        target.completer = completer
        completer.setWidget(target)
        completer.activated.connect(target.insertCompletion)

    def insertCompletion(target, completion):
        tc = target.textCursor()
        extra = (completion.length() -
            target.completer.completionPrefix().length())
        tc.movePosition(QTextCursor.Left)
        tc.movePosition(QTextCursor.EndOfWord)
        tc.insertText(completion.right(extra))
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
            QtCore.Qt.Key_Enter,
            QtCore.Qt.Key_Return,
            QtCore.Qt.Key_Escape,
            QtCore.Qt.Key_Tab,
            QtCore.Qt.Key_Backtab):
                event.ignore()
                return

        ## has ctrl-E been pressed??
        isShortcut = (event.modifiers() == QtCore.Qt.ControlModifier and
                      event.key() == QtCore.Qt.Key_E)
        if (not target.completer or not isShortcut):
            QTextEdit.keyPressEvent(target, event)

        ## ctrl or shift key on it's own??
        ctrlOrShift = event.modifiers() in (QtCore.Qt.ControlModifier ,
                    QtCore.Qt.ShiftModifier)
        if ctrlOrShift and event.text().isEmpty():
            # ctrl or shift key on it's own
            return

        eow = QtCore.QString("~!@#$%^&*()_+{}|:\"<>?,./;'[]\\-=") #end of word

        hasModifier = ((event.modifiers() != QtCore.Qt.NoModifier) and
                        not ctrlOrShift)

        completionPrefix = target.textUnderCursor()

        if (not isShortcut and (hasModifier or event.text().isEmpty() or
        completionPrefix.length() < 2 or
        eow.contains(event.text().right(1)))):
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
