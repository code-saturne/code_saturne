# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2026 EDF S.A.
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
This module defines a push button opening the formula editor, and the
box displaying the current formula next to it.

This module defines the following classes:
- FormulaDisplay
- FormulaButton
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import logging
import math

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

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("FormulaButton")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Formula display
#-------------------------------------------------------------------------------

# The state is shown through a dynamic property and a style sheet rather
# than through the palette, so that it remains visible when a GUI theme
# style sheet is loaded.

_button_style_sheet = """
FormulaButton[formulaState="required"] { color: #d0312d; }
"""

_display_style_sheet = """
FormulaDisplay[formulaState="required"] { color: #d0312d; }
"""

# Width reserved for the formula text in the size hint of a button
# displaying the formula itself, in pixels. Longer formulas are elided,
# the full formula is in the tooltip.

_text_width = 200

# Maximum number of lines shown by the formula box. The box grows with
# the formula up to this number of lines, and scrolls beyond.

_max_lines = 10

# Maximum width of the formula box, in pixels. Longer lines are wrapped.

_max_width = 500


def formula_summary(expression):
    """
    Return a formula on a single line, empty lines and comments
    being removed.
    """
    parts = []
    for line in str(expression).splitlines():
        line = line.split('#', 1)[0].strip()
        if line:
            parts.append(line)
    return " ".join(parts)


def formula_text(expression):
    """
    Return a formula as written by the user, empty lines being removed.
    """
    lines = [l.rstrip() for l in str(expression).splitlines()]
    return "\n".join([l for l in lines if l.strip()])


def _find_layout(layout, widget):
    """
    Return the layout directly containing a widget, searching
    recursively from a given layout.
    """
    if layout is None:
        return None
    if layout.indexOf(widget) >= 0:
        return layout
    for i in range(layout.count()):
        sub = layout.itemAt(i).layout()
        if sub is not None:
            found = _find_layout(sub, widget)
            if found is not None:
                return found
    return None


def _insert_beside(anchor, widget, after):
    """
    Insert a widget just before or after another one in its layout.
    In a grid layout, both widgets share the cell of the anchor.
    Return False if the layout is not handled.
    """
    parent = anchor.parentWidget()
    if parent is None:
        return False
    layout = _find_layout(parent.layout(), anchor)

    if isinstance(layout, QBoxLayout):
        index = layout.indexOf(anchor)
        if after:
            index += 1
        layout.insertWidget(index, widget)
        return True

    if isinstance(layout, QGridLayout):
        index = layout.indexOf(anchor)
        row, column, row_span, column_span = layout.getItemPosition(index)
        item = layout.itemAt(index)
        alignment = item.alignment()
        layout.removeWidget(anchor)
        box = QHBoxLayout()
        box.setContentsMargins(0, 0, 0, 0)
        if after:
            box.addWidget(anchor, 0, alignment)
            box.addWidget(widget)
        else:
            box.addWidget(widget)
            box.addWidget(anchor, 0, alignment)
        layout.addLayout(box, row, column, row_span, column_span)
        return True

    return False


class FormulaDisplay(QTextEdit):
    """
    Read-only box displaying a formula. Its height follows the formula,
    up to a maximum number of lines. Clicking it opens the formula
    editor through the associated button.
    """

    def __init__(self, button, parent=None):
        """
        Constructor.
        """
        QTextEdit.__init__(self, parent)

        self._button = button
        line_edit_hint = QLineEdit().sizeHint()
        self._min_height = line_edit_hint.height()
        self._hint_width = line_edit_hint.width()

        self.setReadOnly(True)
        self.setAcceptRichText(False)
        self.setLineWrapMode(QTextEdit.LineWrapMode.WidgetWidth)
        self.setWordWrapMode(QTextOption.WrapMode.WrapAtWordBoundaryOrAnywhere)
        self.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        self.setTextInteractionFlags(Qt.TextInteractionFlag.NoTextInteraction)
        self.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        self.setSizePolicy(QSizePolicy.Policy.Preferred,
                           QSizePolicy.Policy.Fixed)
        self.viewport().setCursor(Qt.CursorShape.PointingHandCursor)
        self.document().setDocumentMargin(2)

        self.setProperty("formulaState", "none")
        self.setStyleSheet(_display_style_sheet)

        self.__updateHeight()


    def setFormula(self, state, text, tooltip):
        """
        Display a formula with its state ('defined' or 'required').
        """
        if self.property("formulaState") != state:
            self.setProperty("formulaState", state)
            self.style().unpolish(self)
            self.style().polish(self)

        self.setPlainText(text)
        self.setToolTip(tooltip)
        self.updateGeometry()
        self.__updateHeight()


    def __documentHeight(self, width):
        """
        Return the height of the formula wrapped to a given width.
        """
        doc = self.document().clone()
        doc.setTextWidth(width)
        return doc.size().height()


    def __updateHeight(self):
        """
        Adjust the height to the formula, from one line (the height of
        a line edit) up to the maximum number of lines.
        """
        frame = 2*self.frameWidth()
        margins = self.contentsMargins()
        width = self.width() - frame - margins.left() - margins.right()
        height = self.__documentHeight(max(width, 1))

        doc_margin = self.document().documentMargin()
        max_height = _max_lines * self.fontMetrics().lineSpacing() \
                     + 2*doc_margin
        height = min(height, max_height)
        height = int(math.ceil(height)) + frame \
                 + margins.top() + margins.bottom()
        height = max(height, self._min_height)
        if height != self.minimumHeight() or height != self.maximumHeight():
            self.setFixedHeight(height)


    def __textWidth(self):
        """
        Return the width needed to display the longest formula line.
        """
        fm = self.fontMetrics()
        width = 0
        for line in self.toPlainText().splitlines():
            width = max(width, fm.horizontalAdvance(line))
        doc_margin = self.document().documentMargin()
        return int(math.ceil(width + 2*doc_margin)) + 2*self.frameWidth() + 6


    def sizeHint(self):
        """
        Width of the longest formula line, between the width of a line
        edit and a maximum width, longer lines being wrapped.
        """
        hint = QTextEdit.sizeHint(self)
        width = min(max(self.__textWidth(), self._hint_width), _max_width)
        hint.setWidth(width)
        hint.setHeight(self.height())
        return hint


    def minimumSizeHint(self):
        """
        The box may shrink, the formula being wrapped.
        """
        hint = QTextEdit.minimumSizeHint(self)
        hint.setWidth(60)
        hint.setHeight(self._min_height)
        return hint


    def resizeEvent(self, event):
        """
        The formula is wrapped to the new width.
        """
        QTextEdit.resizeEvent(self, event)
        self.__updateHeight()


    def showEvent(self, event):
        """
        The font may change when the box is polished before being shown.
        """
        QTextEdit.showEvent(self, event)
        self.updateGeometry()
        self.__updateHeight()


    def changeEvent(self, event):
        """
        Adjust the size to a new font or style.
        """
        QTextEdit.changeEvent(self, event)
        if event.type() in (QEvent.Type.FontChange, QEvent.Type.StyleChange):
            self.updateGeometry()
            self.__updateHeight()


    def mousePressEvent(self, event):
        """
        Accept the press, the editor being opened on release.
        """
        event.accept()


    def mouseReleaseEvent(self, event):
        """
        Open the formula editor.
        """
        if event.button() == Qt.MouseButton.LeftButton \
           and self._button.isEnabled():
            self._button.click()
        event.accept()


    def wheelEvent(self, event):
        """
        Let the page scroll when the whole formula is displayed.
        """
        if self.verticalScrollBar().maximum() == 0:
            event.ignore()
            return
        QTextEdit.wheelEvent(self, event)

#-------------------------------------------------------------------------------
# Formula button
#-------------------------------------------------------------------------------

class FormulaButton(QPushButton):
    """
    Push button opening the formula editor. The current formula is
    displayed in a box next to the button, which grows with the formula.

    The box is placed just before the button. If the page has a value
    field which is not used when the value is defined by a formula, the
    box may take its place instead (setDisplayField).

    The button has three states:
    - no formula (clearFormula): no box is displayed, as when the value
      is not defined by a formula;
    - formula defined (setFormula): the formula is displayed;
    - formula required (setFormulaRequired): a formula must be defined,
      which is shown in red.

    The constructor only takes a parent, so the class may be used as a
    promoted widget in Qt Designer forms.
    """

    def __init__(self, parent=None):
        """
        Constructor.
        """
        QPushButton.__init__(self, parent)

        self._summary = ""
        self._text = ""
        self._tooltip = ""
        self._default_tooltip = None

        # Formula box, created when a formula is first displayed.
        # If it can not be placed in the layout, the formula is
        # displayed on the button itself.
        self._display = None
        self._on_button = False

        self._field = None
        self._field_hidden = False

        self.setProperty("formulaState", "none")
        self.setStyleSheet(_button_style_sheet)


    def setDisplayField(self, field):
        """
        Display the formula in place of a value field, which is hidden
        while the value is defined by a formula.
        Must be called before the first formula is displayed.
        """
        if self._field is not None:
            self._field.removeEventFilter(self)
        self._field = field
        field.installEventFilter(self)


    def __saveDefaultToolTip(self):
        """
        Keep the tooltip defined in the form. Forms set it after
        construction, so it is saved on first use.
        """
        if self._default_tooltip is None:
            self._default_tooltip = self.toolTip()


    def __display(self):
        """
        Return the formula box, creating it if needed, or None if the
        formula is displayed on the button.
        """
        if self._display is None and not self._on_button:
            anchor = self._field if self._field is not None else self
            display = FormulaDisplay(self, anchor.parentWidget())
            display.hide()
            if _insert_beside(anchor, display, self._field is not None):
                self._display = display
            else:
                display.deleteLater()
                self._on_button = True
                self.updateGeometry()

        return self._display


    def __setState(self, state, summary, text, tooltip):
        """
        Update the state property, the displayed formula and the tooltip.
        """
        if self.property("formulaState") != state:
            self.setProperty("formulaState", state)
            self.style().unpolish(self)
            self.style().polish(self)

        self._summary = summary
        self._text = text
        self._tooltip = tooltip

        display = None
        if state != "none":
            display = self.__display()
        elif self._display is not None:
            display = self._display

        if display is not None:
            if state != "none":
                display.setFormula(state, text, tooltip)
            self.setToolTip(self._default_tooltip)
        elif state != "none":
            self.setToolTip(tooltip)
        else:
            self.setToolTip(self._default_tooltip)

        self.__updateText()
        self.__updateDisplay()


    def setFormula(self, expression):
        """
        Display a defined formula. An empty formula is handled as
        a required one.
        """
        self.__saveDefaultToolTip()

        if not expression:
            self.setFormulaRequired()
            return

        self.__setState("defined",
                        formula_summary(expression),
                        formula_text(expression),
                        str(expression))


    def setFormulaRequired(self, expression=None):
        """
        Mark the formula as required (not defined by the user yet).
        A default expression may be given for display.
        """
        self.__saveDefaultToolTip()

        if expression:
            summary = formula_summary(expression)
            text = formula_text(expression)
            tooltip = str(expression)
        else:
            summary = self.tr("Define formula")
            text = summary
            tooltip = self.tr("Click to define the formula")

        self.__setState("required", summary, text, tooltip)


    def clearFormula(self):
        """
        Display no formula (the value is not defined by a formula).
        """
        self.__saveDefaultToolTip()
        self.__setState("none", "", "", "")


    def formulaState(self):
        """
        Return the current state: 'none', 'defined' or 'required'.
        """
        return self.property("formulaState")


    def formulaDisplay(self):
        """
        Return the box displaying the formula, or None if no formula
        was displayed yet or if it is displayed on the button.
        """
        return self._display


    def __updateDisplay(self):
        """
        Show the formula box when a formula is displayed and the button
        is shown, hiding the associated value field.
        """
        if self._display is None:
            return

        show = self.formulaState() != "none" and not self.isHidden()
        self._display.setVisible(show)
        self._display.setEnabled(self.isEnabled())

        field = self._field
        if field is None:
            return
        try:
            if show and not field.isHidden():
                self._field_hidden = True
                field.hide()
            elif not show and self._field_hidden:
                self._field_hidden = False
                field.show()
        except RuntimeError:
            # The field was deleted with its page
            self._field = None


    def setVisible(self, visible):
        """
        The formula box is shown and hidden with the button.
        """
        QPushButton.setVisible(self, visible)
        self.__updateDisplay()


    def changeEvent(self, event):
        """
        The formula box is enabled and disabled with the button.
        """
        QPushButton.changeEvent(self, event)
        if event.type() == QEvent.Type.EnabledChange:
            self.__updateDisplay()


    def eventFilter(self, obj, event):
        """
        Keep the value field hidden if the page shows it while the
        formula box is displayed in its place.
        """
        if obj is self._field and event.type() == QEvent.Type.Show \
           and self._display is not None and self._display.isVisibleTo(
               self._display.parentWidget()):
            QTimer.singleShot(0, self.__updateDisplay)
        return QPushButton.eventFilter(self, obj, event)


    def __textWidth(self):
        """
        Return the width available for the text.
        """
        w = self.width() - 16
        if not self.icon().isNull():
            w -= self.iconSize().width() + 4
        return max(w, 0)


    def __updateText(self):
        """
        Display the formula summary on the button, elided to the
        available width, when there is no formula box.
        """
        if not self._on_button or not self._summary:
            QPushButton.setText(self, "")
            return

        fm = self.fontMetrics()
        text = fm.elidedText(self._summary,
                             Qt.TextElideMode.ElideRight,
                             self.__textWidth())
        QPushButton.setText(self, text)


    def sizeHint(self):
        """
        Size hint large enough for a short formula when it is displayed
        on the button.
        """
        hint = QPushButton.sizeHint(self)
        if self._on_button:
            w = _text_width + 16
            if not self.icon().isNull():
                w += self.iconSize().width() + 4
            hint.setWidth(max(hint.width(), w))
        return hint


    def minimumSizeHint(self):
        """
        The button may shrink down to its icon.
        """
        hint = QPushButton.minimumSizeHint(self)
        if not self.icon().isNull():
            hint.setWidth(self.iconSize().width() + 16)
        return hint


    def resizeEvent(self, event):
        """
        Elide the formula to the new width.
        """
        QPushButton.resizeEvent(self, event)
        self.__updateText()

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
