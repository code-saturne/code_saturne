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
This module defines the following classes:
- QFileEditor
"""

#-------------------------------------------------------------------------------
# Standard modules
#-------------------------------------------------------------------------------

import sys, os, shutil
from code_saturne import get_cs_data_path
from code_saturne.gui.base import QtGui, QtCore, QtWidgets

from code_saturne.gui.base.QtWidgets import QMainWindow, QMessageBox, \
    QFileDialog, QTextEdit, QPlainTextEdit, QSizePolicy, QMenu, QMessageBox
from code_saturne.gui.base.QtWidgets import QAction

import code_saturne.gui.base.resource_base_rc
from code_saturne.gui.base.SearchBar import SearchBar

#-------------------------------------------------------------------------------
# Optional Pygments import — used for syntax highlighting when available.
# Falls back gracefully to the built-in regex highlighter if not installed.
#-------------------------------------------------------------------------------

try:
    from pygments import lex
    from pygments.lexers import get_lexer_for_filename, guess_lexer, TextLexer
    from pygments.lexers import CLexer, CppLexer, PythonLexer
    from pygments.token import Token
    _PYGMENTS_AVAILABLE = True
except ImportError:
    _PYGMENTS_AVAILABLE = False

#-------------------------------------------------------------------------------
# Compatibility helper: detect which Qt binding is used
#-------------------------------------------------------------------------------

def _fromUtf8(s):
    return s

# Detect Qt binding for compatibility
try:
    from PySide6 import QtCore as _qtc
    _QT_BINDING = 'pyside6'
except ImportError:
    try:
        from PyQt6 import QtCore as _qtc
        _QT_BINDING = 'pyqt6'
    except ImportError:
        _QT_BINDING = 'pyqt5'

def _is_qt6():
    return _QT_BINDING in ('pyside6', 'pyqt6')

#-------------------------------------------------------------------------------
# Local constants
#-------------------------------------------------------------------------------

_tab_size = 2

#-------------------------------------------------------------------------------
# Local functions and/or definitions
#-------------------------------------------------------------------------------

def loc_format(color, style=''):
    """
    Returns a TextCharFormat with the proper attributes
    """

    c = QtGui.QColor()
    c.setNamedColor(color)

    f = QtGui.QTextCharFormat()
    f.setForeground(c)

    # Bold font
    if 'bold' in style:
        f.setFontWeight(QtGui.QFont.Weight.Bold)

    # Italic font
    if 'italic' in style:
        f.setFontItalic(True)

    return f

format_styles = {'keyword'    : loc_format('blue', 'bold'),
                 'operator'   : loc_format('red', 'bold'),
                 'brace'      : loc_format('orange', 'bold'),
                 'string'     : loc_format('magenta', 'italic'),
                 'comment'    : loc_format('darkGreen', 'italic'),
                 'expression' : loc_format('black'),
                 'cs_keyword' : loc_format('darkCyan', 'bold')}

#-------------------------------------------------------------------------------
# Code_Saturne-specific keywords (applied as a second pass in both backends)
#-------------------------------------------------------------------------------

_CS_KEYWORDS = [
    'cs_real_t', 'cs_lnum_t', 'cs_real_3_t',
    'bft_printf', 'bft_printf_flush', 'bft_error',
    'BEGIN_C_DECLS', 'END_C_DECLS',
]

#-------------------------------------------------------------------------------
# HighlightingRule class  (used by the regex fallback backend only)
#-------------------------------------------------------------------------------

class HighlightingRule():

    # ---------------------------------------------------------------
    def __init__(self, pattern, format):

        self.pattern = pattern
        self.format  = format
    # ---------------------------------------------------------------


#-------------------------------------------------------------------------------
# CodeEditor with line numbering
#-------------------------------------------------------------------------------

class LineNumberArea(QtWidgets.QWidget):

    def __init__(self, editor):
        super(LineNumberArea, self).__init__(editor)
        self.editor = editor

    def sizeHint(self):
        return QtCore.QSize(self.editor.lineNumberAreaWidth(),0)

    def paintEvent(self, event):
        self.editor.lineNumberAreaPaintEvent(event)

class CodeEditor(QPlainTextEdit):
    def __init__(self):
        super(CodeEditor, self).__init__()

        self.lineNumberArea = LineNumberArea(self)

        self.blockCountChanged.connect(self.updateLineNumberAreaWidth)
        self.updateRequest.connect(self.updateLineNumberArea)
        self.cursorPositionChanged.connect(self.highlightCurrentLine)

        self.updateLineNumberAreaWidth(0)


    def lineNumberAreaWidth(self):
        digits = 1
        count = max(1, self.blockCount())
        while count >= 10:
            count /= 10
            digits += 1
        # horizontalAdvance replaces width() which is removed in Qt6
        try:
            space = 3 + self.fontMetrics().horizontalAdvance('9') * digits
        except AttributeError:
            space = 3 + self.fontMetrics().width('9') * digits
        return space


    def updateLineNumberAreaWidth(self, _):
        self.setViewportMargins(self.lineNumberAreaWidth(), 0, 0, 0)


    def updateLineNumberArea(self, rect, dy):

        if dy:
            self.lineNumberArea.scroll(0, dy)
        else:
            self.lineNumberArea.update(0, rect.y(), self.lineNumberArea.width(),
                       rect.height())

        if rect.contains(self.viewport().rect()):
            self.updateLineNumberAreaWidth(0)


    def resizeEvent(self, event):
        super(CodeEditor, self).resizeEvent(event)

        cr = self.contentsRect();
        self.lineNumberArea.setGeometry(QtCore.QRect(cr.left(), cr.top(),
                    self.lineNumberAreaWidth(), cr.height()))


    def lineNumberAreaPaintEvent(self, event):
        mypainter = QtGui.QPainter(self.lineNumberArea)

        # Qt.lightGray -> Qt.GlobalColor.lightGray in Qt6
        try:
            mypainter.fillRect(event.rect(), QtCore.Qt.GlobalColor.lightGray)
        except AttributeError:
            mypainter.fillRect(event.rect(), QtCore.Qt.lightGray)

        block = self.firstVisibleBlock()
        blockNumber = block.blockNumber()
        top = self.blockBoundingGeometry(block).translated(self.contentOffset()).top()
        bottom = top + self.blockBoundingRect(block).height()

        x = int(0)
        h = int(self.fontMetrics().height())

        # Qt.AlignRight -> Qt.AlignmentFlag.AlignRight in Qt6
        try:
            flags = int(QtCore.Qt.AlignmentFlag.AlignRight)
        except AttributeError:
            flags = int(QtCore.Qt.AlignRight)

        while block.isValid() and (top <= event.rect().bottom()):
            if block.isVisible() and (bottom >= event.rect().top()):
                number = str(blockNumber + 1)
                # Qt.black -> Qt.GlobalColor.black in Qt6
                try:
                    mypainter.setPen(QtCore.Qt.GlobalColor.black)
                except AttributeError:
                    mypainter.setPen(QtCore.Qt.black)
                y = int(top)
                w = int(self.lineNumberArea.width())
                mypainter.drawText(x, y, w, h, flags, number)

            block = block.next()
            top = bottom
            bottom = top + self.blockBoundingRect(block).height()
            blockNumber += 1


    def highlightCurrentLine(self):
        extraSelections = []

        if not self.isReadOnly():
            selection = QTextEdit.ExtraSelection()

            # Qt.yellow -> Qt.GlobalColor.yellow in Qt6
            try:
                lineColor = QtGui.QColor(QtCore.Qt.GlobalColor.yellow).lighter(160)
            except AttributeError:
                lineColor = QtGui.QColor(QtCore.Qt.yellow).lighter(160)

            selection.format.setBackground(lineColor)
            # QTextFormat.FullWidthSelection -> QTextFormat.Property.FullWidthSelection in Qt6
            try:
                selection.format.setProperty(QtGui.QTextFormat.Property.FullWidthSelection, True)
            except AttributeError:
                selection.format.setProperty(QtGui.QTextFormat.FullWidthSelection, True)
            selection.cursor = self.textCursor()
            selection.cursor.clearSelection()
            extraSelections.append(selection)
        self.setExtraSelections(extraSelections)

#-------------------------------------------------------------------------------
# QtextHighlighter — public facade
#
# Instantiates the Pygments backend when available, the regex backend otherwise.
# The caller (QFileEditor.newFile) does not need to change.
#-------------------------------------------------------------------------------

def QtextHighlighter(parent, extension):
    """
    Factory function returning the best available highlighter for *extension*.

    Returns a QSyntaxHighlighter subclass instance attached to *parent*
    (a QTextDocument). Uses Pygments when installed, falls back to the
    built-in regex highlighter otherwise.
    """
    if _PYGMENTS_AVAILABLE:
        return _PygmentsHighlighter(parent, extension)
    else:
        return _RegexHighlighter(parent, extension)


#-------------------------------------------------------------------------------
# _PygmentsHighlighter — Pygments-based backend
#-------------------------------------------------------------------------------

# Map file extensions to Pygments lexer classes.
# Plain text (unknown extension) gets the null TextLexer so Pygments still
# runs and the CS keyword pass still fires.
_EXT_TO_LEXER = {
    'c'   : None,  # resolved dynamically via get_lexer_for_filename
    'cpp' : None,
    'cxx' : None,
    'c++' : None,
    'h'   : None,
    'hxx' : None,
    'f90' : None,
    'F90' : None,
    'F'   : None,
    'f77' : None,
    'py'  : None,
}

# Map Pygments token types to our format_styles keys.
# Order matters: more specific entries first.
_TOKEN_FORMAT = [
    (Token.Comment,               'comment'),
    (Token.Literal.String,        'string'),
    (Token.Literal.String.Doc,    'comment'),   # docstrings → green italic
    (Token.Keyword,               'keyword'),
    (Token.Keyword.Type,          'keyword'),
    (Token.Name.Builtin,          'keyword'),
    (Token.Operator,              'operator'),
    (Token.Punctuation,           'brace'),
    (Token.Literal.Number,        'expression'),
]


class _PygmentsHighlighter(QtGui.QSyntaxHighlighter):
    """
    Syntax highlighter backed by Pygments.

    For each document block (line), Pygments tokenizes the *entire* document
    text once (cached), then we apply Qt formats token by token for the
    current line.  A second pass applies the CS-specific keywords on top.
    """

    # ---------------------------------------------------------------
    def __init__(self, parent_doc, extension):
        QtGui.QSyntaxHighlighter.__init__(self, parent_doc)
        self.parent_doc = parent_doc
        self.extension  = extension

        # Resolve the Pygments lexer for this extension.
        self._lexer = self._get_lexer(extension)

        # Pre-build CS keyword patterns (QRegularExpression or QRegExp).
        self._cs_patterns = self._build_cs_patterns()
    # ---------------------------------------------------------------

    # ---------------------------------------------------------------
    def _get_lexer(self, extension):
        """Return the best Pygments lexer for *extension*."""
        fake_filename = 'file.' + extension
        try:
            lexer = get_lexer_for_filename(fake_filename, stripnl=False)
        except Exception:
            lexer = TextLexer(stripnl=False)
        return lexer
    # ---------------------------------------------------------------

    # ---------------------------------------------------------------
    def _build_cs_patterns(self):
        """
        Pre-compile regex patterns for Code_Saturne keywords.
        Returns a list of compiled pattern objects (QRegularExpression or QRegExp).
        """
        patterns = []
        for kw in _CS_KEYWORDS:
            if _is_qt6():
                patterns.append(QtCore.QRegularExpression(r'\b' + kw + r'\b'))
            else:
                patterns.append(QtCore.QRegExp(r'\b' + kw + r'\b'))
        return patterns
    # ---------------------------------------------------------------

    # ---------------------------------------------------------------
    def _get_tokens_for_block(self, block_text, block_start):
        """
        Run Pygments on the full document text and extract only the tokens
        that fall inside [block_start, block_start + len(block_text)].

        Returns list of (relative_start, length, format_key) tuples.
        """
        # Retrieve full document text each time.
        # QSyntaxHighlighter gives us no cheap way to cache across blocks
        # without subclassing QTextDocument, so we accept the overhead —
        # acceptable for files up to a few thousand lines.
        doc = self.document()
        full_text = doc.toPlainText()

        results = []
        pos = 0
        block_end = block_start + len(block_text)

        for ttype, value in lex(full_text, self._lexer):
            token_start = pos
            token_end   = pos + len(value)
            pos = token_end

            # Skip tokens entirely outside this block.
            if token_end <= block_start:
                continue
            if token_start >= block_end:
                break

            # Clamp to block boundaries.
            rel_start  = max(token_start, block_start) - block_start
            rel_end    = min(token_end,   block_end)   - block_start
            rel_length = rel_end - rel_start
            if rel_length <= 0:
                continue

            fmt_key = self._resolve_format(ttype)
            if fmt_key:
                results.append((rel_start, rel_length, fmt_key))

        return results
    # ---------------------------------------------------------------

    # ---------------------------------------------------------------
    @staticmethod
    def _resolve_format(ttype):
        """
        Map a Pygments token type to a format_styles key.
        Walks the token hierarchy from specific to generic.
        """
        for base_type, fmt_key in _TOKEN_FORMAT:
            if ttype is base_type or ttype in base_type:
                return fmt_key
        return None
    # ---------------------------------------------------------------

    # ---------------------------------------------------------------
    def highlightBlock(self, text):
        """Apply Pygments-based highlighting, then the CS keyword pass."""

        # Compute absolute offset of this block in the document.
        block       = self.currentBlock()
        block_start = block.position()

        # --- Pygments pass ---
        for rel_start, length, fmt_key in self._get_tokens_for_block(text, block_start):
            self.setFormat(rel_start, length, format_styles[fmt_key])

        # --- Code_Saturne keyword pass (second pass, always overwrites) ---
        self._apply_cs_keywords(text)
    # ---------------------------------------------------------------

    # ---------------------------------------------------------------
    def _apply_cs_keywords(self, text):
        """Highlight Code_Saturne-specific identifiers in *text*."""
        fmt = format_styles['cs_keyword']
        for pat in self._cs_patterns:
            if _is_qt6():
                it = pat.globalMatch(text)
                while it.hasNext():
                    m = it.next()
                    self.setFormat(m.capturedStart(), m.capturedLength(), fmt)
            else:
                idx = pat.indexIn(text)
                while idx >= 0:
                    self.setFormat(idx, pat.matchedLength(), fmt)
                    idx = pat.indexIn(text, idx + pat.matchedLength())
    # ---------------------------------------------------------------


#-------------------------------------------------------------------------------
# _RegexHighlighter — pure-Qt regex fallback backend
#
# Functionally equivalent to the original QtextHighlighter, with:
#   - cleaner rule construction (one helper instead of duplicated Qt5/Qt6 blocks)
#   - Code_Saturne keywords integrated into the C/C++ keyword list
#   - CS keyword second pass (same as the Pygments backend, for consistency)
#-------------------------------------------------------------------------------

class _RegexHighlighter(QtGui.QSyntaxHighlighter):
    """
    Syntax highlighter using QRegularExpression (Qt6) or QRegExp (Qt5).
    Zero external dependencies.
    """

    # ---------------------------------------------------------------
    def __init__(self, parent_doc, extension):
        QtGui.QSyntaxHighlighter.__init__(self, parent_doc)
        self.parent_doc = parent_doc
        self.highlightingRules = []

        # --- Keyword lists ---

        fortran_kw = ['if', 'else', 'endif', 'do', 'enddo', 'end',
                      'implicit none', 'use', 'subroutine', 'function',
                      'double precision', 'real', 'integer', 'char',
                      'allocatable', 'allocate', 'deallocate', 'dimension',
                      'select case', 'call']

        # Standard C/C++ keywords + Code_Saturne-specific identifiers.
        c_kw = ['if', 'else', 'for', 'switch', 'while',
                '\#', 'include', 'pass', 'return', 'del', 'delete',
                'assert', 'true', 'false', 'continue', 'break',
                'fprintf', 'bft_printf', 'bft_printf_flush', 'bft_error',
                'cs_real_t', 'cs_lnum_t', 'cs_real_3_t',
                'int', 'char', 'string', 'void', 'double', 'const',
                'BEGIN_C_DECLS', 'END_C_DECLS']

        py_kw = ['if', 'elif', 'for', 'range', 'while', 'return', 'def',
                 'True', 'False']

        self.kw = []
        if extension in ['f90', 'F90', 'F', 'f77']:
            for kw in fortran_kw:
                self.kw.append(kw)
                self.kw.append(kw.upper())
        elif extension in ['c', 'cpp', 'cxx', 'c++']:
            for kw in c_kw:
                self.kw.append(kw)
                self.kw.append(kw.upper())
        elif extension == 'py':
            for kw in py_kw:
                self.kw.append(kw)

        # Operators
        self.op = ['=', '==', '!=', '<', '>', '<=', '>=',
                   '\+', '-', '\*', '/', '\%', '\*\*',
                   '\+=', '-=', '\*=', '/=', '->', '=>',
                   '\^', '\|', '\&', '\|\|', '\&\&']

        # Braces
        self.br = ['\(', '\)', '\{', '\}', '\[', '\]']

        # Build rules using the appropriate regex class.
        self._build_rules()

        # CS keyword patterns (second pass, same as Pygments backend).
        self._cs_patterns = self._build_cs_patterns()
    # ---------------------------------------------------------------

    # ---------------------------------------------------------------
    def _make_pattern(self, expr):
        """Return a QRegularExpression (Qt6) or QRegExp (Qt5) for *expr*."""
        if _is_qt6():
            return QtCore.QRegularExpression(expr)
        else:
            return QtCore.QRegExp(expr)
    # ---------------------------------------------------------------

    # ---------------------------------------------------------------
    def _build_rules(self):
        """Construct all highlighting rules."""
        R = HighlightingRule

        for kw in self.kw:
            self.highlightingRules.append(
                R(self._make_pattern(r'\b' + kw + r'\b'), format_styles['keyword'])
            )

        for op in self.op:
            self.highlightingRules.append(
                R(self._make_pattern(op), format_styles['operator'])
            )

        for br in self.br:
            self.highlightingRules.append(
                R(self._make_pattern(br), format_styles['brace'])
            )

        # String literals — pattern differs slightly between Qt5 and Qt6
        # due to backslash handling in QRegExp vs QRegularExpression.
        if _is_qt6():
            str_pat  = r'"[^"\\]*(\\.[^"\\]*)*"'
            cmt_pat  = r'//[^\n]*'
            fcmt_pat = r'![^\n]*'
            num_pats = [r'[+-]?[0-9]+[lL]?',
                        r'[+-]?0[xX][0-9A-Fa-f]+[lL]?',
                        r'[+-]?[0-9]+(?:\.[0-9]+)?(?:[eE][+-]?[0-9]+)?']
        else:
            str_pat  = r'"[^"\\]*(\\.[^"\\]*)*"'
            cmt_pat  = r'//[^\n]*'
            fcmt_pat = r'![^\n]*'
            num_pats = [r'[+-]?[0-9]+[lL]?',
                        r'[+-]?0[xX][0-9A-Fa-f]+[lL]?',
                        r'[+-]?[0-9]+(?:\.[0-9]+)?(?:[eE][+-]?[0-9]+)?']

        self.highlightingRules.append(
            R(self._make_pattern(str_pat), format_styles['string'])
        )
        self.highlightingRules.append(
            R(self._make_pattern(cmt_pat), format_styles['comment'])
        )
        self.highlightingRules.append(
            R(self._make_pattern(fcmt_pat), format_styles['comment'])
        )
        for np in num_pats:
            self.highlightingRules.append(
                R(self._make_pattern(np), format_styles['expression'])
            )
    # ---------------------------------------------------------------

    # ---------------------------------------------------------------
    def _build_cs_patterns(self):
        """Pre-compile patterns for the CS keyword second pass."""
        patterns = []
        for kw in _CS_KEYWORDS:
            patterns.append(self._make_pattern(r'\b' + kw + r'\b'))
        return patterns
    # ---------------------------------------------------------------

    # ---------------------------------------------------------------
    def highlightBlock(self, text):
        """Apply regex rules, then the CS keyword pass."""

        for rule in self.highlightingRules:
            if _is_qt6():
                it = rule.pattern.globalMatch(text)
                while it.hasNext():
                    m = it.next()
                    index  = m.capturedStart()
                    length = m.capturedLength()
                    self.setFormat(index, length, rule.format)
            else:
                exp   = QtCore.QRegExp(rule.pattern)
                index = exp.indexIn(text)

                while index >= 0:
                    length = exp.matchedLength()
                    ok_to_highlight = True
                    if len(text) > index + length:
                        if text[index + length] not in self.op + [' ']:
                            ok_to_highlight = False
                    if text[index:index + length] not in self.op + self.br:
                        ok_to_highlight = True

                    if ok_to_highlight:
                        self.setFormat(index, length, rule.format)
                    index = text.find(exp.cap(), index + length)

        self.setCurrentBlockState(0)

        # C/C++ multi-line comments  /* … */
        self.highlightCommentsOverLines(text, r'/\*', r'\*/')

        # CS keyword second pass.
        self._apply_cs_keywords(text)
    # ---------------------------------------------------------------

    # ---------------------------------------------------------------
    def _apply_cs_keywords(self, text):
        """Highlight Code_Saturne-specific identifiers (second pass)."""
        fmt = format_styles['cs_keyword']
        for pat in self._cs_patterns:
            if _is_qt6():
                it = pat.globalMatch(text)
                while it.hasNext():
                    m = it.next()
                    self.setFormat(m.capturedStart(), m.capturedLength(), fmt)
            else:
                idx = pat.indexIn(text)
                while idx >= 0:
                    self.setFormat(idx, pat.matchedLength(), fmt)
                    idx = pat.indexIn(text, idx + pat.matchedLength())
    # ---------------------------------------------------------------

    # ---------------------------------------------------------------
    def highlightCommentsOverLines(self, text, dls, dle):
        """Handle C/C++ block comments spanning multiple lines."""

        if _is_qt6():
            startExpression = QtCore.QRegularExpression(dls)
            endExpression   = QtCore.QRegularExpression(dle)
            ref_state = 1

            if self.previousBlockState() == ref_state:
                start = 0
                add   = 0
            else:
                m     = startExpression.match(text)
                start = m.capturedStart() if m.hasMatch() else -1
                add   = m.capturedLength() if m.hasMatch() else 0

            while start >= 0:
                m_end = endExpression.match(text, start + add)
                end   = m_end.capturedStart() if m_end.hasMatch() else -1

                if end >= add:
                    length = end - start + add + m_end.capturedLength()
                    self.setCurrentBlockState(0)
                else:
                    self.setCurrentBlockState(ref_state)
                    length = len(text) - start + add

                self.setFormat(start, length, format_styles['comment'])
                m_next = endExpression.match(text, start + length)
                start  = m_next.capturedStart() if m_next.hasMatch() else -1
        else:
            startExpression = QtCore.QRegExp(dls)
            endExpression   = QtCore.QRegExp(dle)
            ref_state = 1

            if self.previousBlockState() == ref_state:
                start = 0
                add   = 0
            else:
                start = startExpression.indexIn(text)
                add   = startExpression.matchedLength()

            while start >= 0:
                end = endExpression.indexIn(text, start + add)

                if end >= add:
                    length = end - start + add + endExpression.matchedLength()
                    self.setCurrentBlockState(0)
                else:
                    self.setCurrentBlockState(ref_state)
                    length = len(text) - start + add

                self.setFormat(start, length, format_styles['comment'])
                start = endExpression.indexIn(text, start + length)
    # ---------------------------------------------------------------


#-------------------------------------------------------------------------------
# QMessageBox which expands
#-------------------------------------------------------------------------------

class QExpandingMessageBox(QMessageBox):
    """
    A QMessageBox which expands.
    """

    def __init__(self, parent=None):
        QMessageBox.__init__(self,parent=parent)
        self.setSizeGripEnabled(True)

    def event(self, ev):

        result = QMessageBox.event(self, ev)

        self.setMinimumHeight(10)
        self.setMaximumHeight(16777215)
        self.setMinimumWidth(10)
        self.setMaximumWidth(16777215)
        self.setSizePolicy(QSizePolicy.Policy.Expanding,
                           QSizePolicy.Policy.Expanding)

        text = self.findChild(QTextEdit)
        if text != None:
            self.setMinimumHeight(10)
            self.setMaximumHeight(16777215)
            self.setMinimumWidth(1050)
            self.setMaximumWidth(16777215)

            text.setMinimumHeight(10)
            text.setMaximumHeight(16777215)
            text.setMinimumWidth(1000)
            text.setMaximumWidth(16777215)
            text.setSizePolicy(QSizePolicy.Policy.Expanding,
                               QSizePolicy.Policy.Expanding)

        return result

#-------------------------------------------------------------------------------
# QFileSystemModel with modified header
#-------------------------------------------------------------------------------

class FileSystemModel(QtWidgets.QFileSystemModel):

    def __init__(self, title):
        """
        """
        QtWidgets.QFileSystemModel.__init__(self)
        self.title = title


    def headerData(self, section, orientation, role):
        if orientation == QtCore.Qt.Orientation.Horizontal \
          and role == QtCore.Qt.ItemDataRole.DisplayRole:
            if section == 0:
                return self.tr(self.title)
        return None


#-------------------------------------------------------------------------------
# Helper for QMessageBox Yes/No compatibility
#-------------------------------------------------------------------------------

def _msgbox_yes():
    try:
        return QMessageBox.StandardButton.Yes
    except AttributeError:
        return QMessageBox.Yes

def _msgbox_no():
    try:
        return QMessageBox.StandardButton.No
    except AttributeError:
        return QMessageBox.No

def _msgbox_yes_no():
    return _msgbox_yes() | _msgbox_no()


#-------------------------------------------------------------------------------
# Explorer class
#-------------------------------------------------------------------------------

class Explorer():
    """
    Editor class. Used for file editing and/or viewing
    """

    # ---------------------------------------------------------------
    def __init__(self, parent=None, root_dir=None, dir_type=None,
                 case_name=None, readOnly=False):

        self.parent_widget = parent  # renamed to avoid collision with Qt parent()

        self.root_dir = root_dir
        self.dir_type = dir_type

        self.readOnly = readOnly
        self.readerMode = readOnly

        # Explorer
        self.explorer = self._initFileExplorer()
        self._initExplorerActions(case_name)

    # ---------------------------------------------------------------


    # ---------------------------------------------------------------
    def _initFileExplorer(self):
        """
        Create the File explorer object based on the QFileSystemModel widget.
        """

        if self.dir_type == 'SHARE':
            name = 'Reference'
        elif self.dir_type in ('SRC', 'DATA'):
            name = 'User files'
        else:
            name = 'Name'

        model = FileSystemModel(name)
        if self.root_dir:
            model.setRootPath(self.root_dir)

        tree = QtWidgets.QTreeView(None)

        tree.setModel(model)
        tree.setSortingEnabled(True)
        tree.setWindowTitle('Explorer')
        if self.root_dir:
            tree.setRootIndex(model.index(self.root_dir))

        # Hide unnecessary columns
        nc = tree.header().count()

        for i in range(1, nc):
            tree.hideColumn(i)

        # Right click menu
        tree.setContextMenuPolicy(QtCore.Qt.ContextMenuPolicy.CustomContextMenu)
        tree.customContextMenuRequested.connect(self.explorerContextMenu)

        # Double click
        tree.doubleClicked.connect(self._explorerDoubleClick)

        return tree;
    # ---------------------------------------------------------------


    # ---------------------------------------------------------------
    def _initExplorerActions(self, case_name=None):
        """
        Create explorer actions dictionary
        """

        if case_name:
            case_dir_name = str(case_name)
        else:
            case_dir_name = 'SRC'

        _editAction = QAction(self.explorer.model())
        _editAction.setText('Edit file')
        _editAction.triggered.connect(self.parent_widget._editSelectedFile)

        _viewAction = QAction(self.explorer.model())
        _viewAction.setText('View file')
        _viewAction.triggered.connect(self.parent_widget._viewSelectedFile)

        _copyAction = QAction(self.explorer.model())
        _copyAction.setText('Copy to ' + case_dir_name)
        _copyAction.triggered.connect(self.parent_widget._copySelectedFile)

        _removeAction = QAction(self.explorer.model())
        _removeAction.setText('Remove from ' + case_dir_name)
        _removeAction.triggered.connect(self.parent_widget._removeSelectedFile)

        _restoreAction = QAction(self.explorer.model())
        _restoreAction.setText('Move to ' + case_dir_name)
        _restoreAction.triggered.connect(self.parent_widget._restoreSelectedFile)

        _deleteAction = QAction(self.explorer.model())
        _deleteAction.setText('Delete')
        _deleteAction.triggered.connect(self.parent_widget._deleteSelectedFile)

        self._explorerActions = {'edit': _editAction,
                                 'view': _viewAction,
                                 'copy': _copyAction,
                                 'remove': _removeAction,
                                 'restore': _restoreAction,
                                 'delete': _deleteAction}
    # ---------------------------------------------------------------


    # ---------------------------------------------------------------
    def _updateCurrentSelection(self):
        """
        Update the current selection
        """
        path2file = ''
        fname = ''
        for idx in self.explorer.selectedIndexes():
            fname = idx.data(QtCore.Qt.ItemDataRole.DisplayRole)
            c = idx
            p = c.parent()
            ps = p.data(QtCore.Qt.ItemDataRole.DisplayRole)
            while True:
                ctxt = c.data(QtCore.Qt.ItemDataRole.DisplayRole)
                ptxt = p.data(QtCore.Qt.ItemDataRole.DisplayRole)
                if ptxt in [None, self.parent_widget.case_name]:
                    pe = ptxt
                    break
                path2file = ptxt + '/' + path2file
                c = p
                p = c.parent()

        if len(fname) > 0:
            self.parent_widget._currentSelection = {'filename':fname,
                                             'subpath' :path2file,
                                             'filedir' :ps,
                                             'origdir' :pe}

        return
    # ---------------------------------------------------------------


    # ---------------------------------------------------------------
    def _explorerDoubleClick(self):
        """
        Double click action
        """

        self._updateCurrentSelection()

        clicked = os.path.join(self.parent_widget._currentSelection['subpath'],
                               self.parent_widget._currentSelection['filename'])

        if self.root_dir:
            clicked = os.path.join(self.root_dir, clicked)

        edit_list = ['SRC']

        if not os.path.isdir(clicked):
            if self.parent_widget._currentSelection['filedir'] in edit_list:
                self.parent_widget._editSelectedFile()
            else:
                self.parent_widget._viewSelectedFile()

    # ---------------------------------------------------------------


    # ---------------------------------------------------------------
    def explorerContextMenu(self, position):
        """
        Custom menu for the mouse right-click.
        """

        self._updateCurrentSelection()

        path2file = self.parent_widget._currentSelection['subpath']
        fname     = self.parent_widget._currentSelection['filename']
        pe        = self.parent_widget._currentSelection['origdir']
        ps        = self.parent_widget._currentSelection['filedir']

        self._contextMenu = QMenu()

        if (path2file == '' or path2file is None ) and self.root_dir:
            path2file = self.root_dir

        if self.dir_type == 'SHARE':
            if not os.path.isdir(os.path.join(path2file, fname)):
                self._contextMenu.addAction(self._explorerActions['view'])
                self._contextMenu.addAction(self._explorerActions['copy'])
        elif pe == 'SRC':
            if not os.path.isdir(os.path.join(path2file, fname)):
                if ps == 'SRC':
                    self._contextMenu.addAction(self._explorerActions['edit'])
                    self._contextMenu.addAction(self._explorerActions['remove'])
                elif ps in ['EXAMPLES', 'REFERENCE']:
                    self._contextMenu.addAction(self._explorerActions['view'])
                    self._contextMenu.addAction(self._explorerActions['copy'])
                elif ps in ['DRAFT']:
                    self._contextMenu.addAction(self._explorerActions['view'])
                    self._contextMenu.addAction(self._explorerActions['restore'])
                    self._contextMenu.addAction(self._explorerActions['delete'])
        elif pe == 'DATA':
            if not os.path.isdir(os.path.join(path2file, fname)):
                if ps == 'DATA':
                    if fname not in ('setup.xml', 'run.cfg'):
                        self._contextMenu.addAction(self._explorerActions['edit'])
                        self._contextMenu.addAction(self._explorerActions['remove'])
                    else:
                        self._contextMenu.addAction(self._explorerActions['view'])
                elif ps in ['REFERENCE']:
                    self._contextMenu.addAction(self._explorerActions['view'])
                    self._contextMenu.addAction(self._explorerActions['copy'])
                elif ps in ['DRAFT']:
                    self._contextMenu.addAction(self._explorerActions['view'])
                    self._contextMenu.addAction(self._explorerActions['restore'])
                    self._contextMenu.addAction(self._explorerActions['delete'])
        else:
            if not os.path.isdir(os.path.join(path2file, fname)):
                self._contextMenu.addAction(self._explorerActions['view'])

        # exec_() -> exec() in Qt6
        try:
            self._contextMenu.exec(self.explorer.viewport().mapToGlobal(position))
        except AttributeError:
            self._contextMenu.exec_(self.explorer.viewport().mapToGlobal(position))
    # ---------------------------------------------------------------


#-------------------------------------------------------------------------------
# QFileEditor class
#-------------------------------------------------------------------------------

class QFileEditor(QMainWindow):
    """
    Editor class. Used for file editing and/or viewing
    """

    # ---------------------------------------------------------------
    def __init__(self, parent=None, case_dir=None, reference_dir=None,
                 readOnly=False, noOpen=False, useHighlight=True):
        super(QFileEditor, self).__init__(parent)
        self.setGeometry(50, 50, 500, 300)

        self.setWindowTitle("code_saturne built-in file editor")
        self.parent_widget = parent  # renamed to avoid collision with Qt parent()

        self.case_dir = case_dir
        if self.case_dir:
            self.case_name = os.path.split(case_dir)[-1]

        self.last_dir = case_dir

        self.readOnly = readOnly
        self.readerMode = readOnly

        self.useHighlight = useHighlight

        self.opened = False
        self.saved  = True

        # Path to icons:
        icons_path = os.path.join(get_cs_data_path(), 'icons', '22x22')

        # Open file action
        open_img_path = os.path.join(icons_path, 'document-open.png')
        icon_open     = QtGui.QIcon()
        icon_open.addPixmap(QtGui.QPixmap(_fromUtf8(open_img_path)),
                            QtGui.QIcon.Mode.Normal,
                            QtGui.QIcon.State.Off)
        self.openFileAction = QAction(icon_open, "Open", self)
        self.openFileAction.setShortcut("Ctrl+O")
        self.openFileAction.setStatusTip('Open File')
        self.openFileAction.triggered.connect(self.openFileForAction)

        # New file action
        new_img_path = os.path.join(icons_path, 'document-new.png')
        icon_new     = QtGui.QIcon()
        icon_new.addPixmap(QtGui.QPixmap(new_img_path),
                          QtGui.QIcon.Mode.Normal,
                          QtGui.QIcon.State.Off)
        self.newFileAction = QAction(icon_new, "New", self)
        self.newFileAction.setShortcut("Ctrl+E")
        self.newFileAction.setStatusTip('Create new file')
        self.newFileAction.triggered.connect(self.newFile)

        # Save action
        save_img_path = os.path.join(icons_path, 'document-save.png')
        icon_save     = QtGui.QIcon()
        icon_save.addPixmap(QtGui.QPixmap(save_img_path),
                          QtGui.QIcon.Mode.Normal,
                          QtGui.QIcon.State.Off)
        self.saveFileAction = QAction(icon_save, "Save", self)
        self.saveFileAction.setShortcut("Ctrl+S")
        self.saveFileAction.setStatusTip('Save file')
        self.saveFileAction.triggered.connect(self.saveFile)

        # Save as action
        saveas_img_path = os.path.join(icons_path, 'document-save-as.png')
        icon_saveas     = QtGui.QIcon()
        icon_saveas.addPixmap(QtGui.QPixmap(saveas_img_path),
                              QtGui.QIcon.Mode.Normal,
                              QtGui.QIcon.State.Off)
        self.saveFileAsAction = QAction(icon_saveas, "Save as", self)
        self.saveFileAsAction.setStatusTip('Save file as')
        self.saveFileAsAction.triggered.connect(self.saveFileAs)

        # Close file action
        close_img_path = os.path.join(icons_path, 'process-stop.png')
        icon_close     = QtGui.QIcon()
        icon_close.addPixmap(QtGui.QPixmap(close_img_path),
                             QtGui.QIcon.Mode.Normal,
                             QtGui.QIcon.State.Off)
        self.closeFileAction = QAction(icon_close, "Close file", self)
        self.closeFileAction.setShortcut("Ctrl+Q")
        self.closeFileAction.setStatusTip('Close opened file')
        self.closeFileAction.triggered.connect(self.closeOpenedFile)

        # Exit editor action
        quit_img_path = os.path.join(icons_path, 'system-log-out.png')
        icon_quit     = QtGui.QIcon()
        icon_quit.addPixmap(QtGui.QPixmap(quit_img_path),
                          QtGui.QIcon.Mode.Normal,
                          QtGui.QIcon.State.Off)
        self.quitAction = QAction(icon_quit, "Quit", self)
        self.quitAction.setStatusTip('Quit the editor')
        self.quitAction.triggered.connect(self.closeApplication)

        self.statusBar()

        # File toolbar
        self.toolbar = self.addToolBar("Options")

        self.toolbar.addAction(self.newFileAction)
        if not noOpen:
            self.toolbar.addAction(self.openFileAction)
        self.toolbar.addAction(self.saveFileAction)
        self.toolbar.addAction(self.saveFileAsAction)
        self.toolbar.addAction(self.closeFileAction)
        self.toolbar.addAction(self.quitAction)

        # File menu
        self.mainMenu = self.menuBar()

        self.fileMenu = self.mainMenu.addMenu('&File')
        self.fileMenu.addAction(self.newFileAction)
        if not noOpen:
            self.fileMenu.addAction(self.openFileAction)
        self.fileMenu.addAction(self.saveFileAction)
        self.fileMenu.addAction(self.saveFileAsAction)
        self.fileMenu.addAction(self.closeFileAction)
        self.fileMenu.addAction(self.quitAction)

        # Explorer
        self.explorer = Explorer(parent=self,
                                 root_dir=self.case_dir,
                                 dir_type=self.case_name,
                                 case_name=self.case_name)

        # DockWidget flags compatibility
        try:
            _no_features = QtWidgets.QDockWidget.DockWidgetFeature.NoDockWidgetFeatures
            _left_area   = QtCore.Qt.DockWidgetArea.LeftDockWidgetArea
            _right_area  = QtCore.Qt.DockWidgetArea.RightDockWidgetArea
            _top_area    = QtCore.Qt.DockWidgetArea.TopDockWidgetArea
            _bottom_area = QtCore.Qt.DockWidgetArea.BottomDockWidgetArea
        except AttributeError:
            _no_features = QtWidgets.QDockWidget.NoDockWidgetFeatures
            _left_area   = QtCore.Qt.LeftDockWidgetArea
            _right_area  = QtCore.Qt.RightDockWidgetArea
            _top_area    = QtCore.Qt.TopDockWidgetArea
            _bottom_area = QtCore.Qt.BottomDockWidgetArea

        dock = QtWidgets.QDockWidget("User files explorer", self)
        dock.setAllowedAreas(_left_area | _right_area)
        dock.setFeatures(_no_features)
        dock.setWidget(self.explorer.explorer)
        self.addDockWidget(_left_area, dock)

        # Explorer ref
        self.explorer_ref = None
        if reference_dir:
            self.explorer_ref = Explorer(parent=self,
                                         root_dir=reference_dir,
                                         dir_type='SHARE',
                                         case_name=self.case_name)
            dock = QtWidgets.QDockWidget("Reference files explorer", self)
            dock.setAllowedAreas(_left_area | _right_area)
            dock.setFeatures(_no_features)
            dock.setWidget(self.explorer_ref.explorer)
            self.addDockWidget(_left_area, dock)

        # Editor
        self.textEdit = self._initFileEditor()
        self.setCentralWidget(self.textEdit)

        # Settings
        settings = QtCore.QSettings()

        try:
            self.restoreGeometry(settings.value("MainWindow/Geometry", QtCore.QByteArray()))
            self.restoreState(settings.value("MainWindow/State", QtCore.QByteArray()))
        except:
            self.recentFiles = settings.value("RecentFiles").toStringList()
            self.restoreGeometry(settings.value("MainWindow/Geometry").toByteArray())
            self.restoreState(settings.value("MainWindow/State").toByteArray())

        # file attributes
        self.filename = ""
        self.file_extension  = ""

        # searchbar
        self.searchBar = SearchBar(self.textEdit)

        dock = QtWidgets.QDockWidget("Editor search bar", self)
        dock.setAllowedAreas(_top_area | _bottom_area)
        dock.setFeatures(_no_features)
        dock.setWidget(self.searchBar)
        self.addDockWidget(_top_area, dock)


    # ---------------------------------------------------------------
    def _initFileEditor(self):
        """
        Create the Editor widget based on QTextEdit
        """

        base_font = QtGui.QFont()
        base_font.setFamily("Courier")
        base_font.setStyleHint(QtGui.QFont.StyleHint.Monospace)
        base_font.setFixedPitch(True)
        base_font.setPointSize(10)

        font_metrics = QtGui.QFontMetrics(base_font)
        _tab_string = ''
        for i in range(_tab_size):
            _tab_string += ' '

        textEdit = CodeEditor()
        textEdit.setFont(base_font)
        textEdit.textChanged.connect(self.updateFileState)
        textEdit.setReadOnly(self.readOnly)
        policy = textEdit.sizePolicy()
        policy.setHorizontalPolicy(QSizePolicy.Policy.Expanding)
        textEdit.setSizePolicy(policy)

        # setTabStopWidth removed in Qt6, use setTabStopDistance
        try:
            textEdit.setTabStopDistance(font_metrics.horizontalAdvance(_tab_string))
        except AttributeError:
            try:
                textEdit.setTabStopWidth(font_metrics.horizontalAdvance(_tab_string))
            except AttributeError:
                textEdit.setTabStopWidth(font_metrics.width(_tab_string))

        return textEdit
    # ---------------------------------------------------------------


    # ---------------------------------------------------------------
    def _initFileExplorer(self, base_dir=None, name="User Files"):
        """
        Create the File explorer object based on the QFileSystemModel widget.
        """

        model = FileSystemModel(name)
        rootp = ''
        if base_dir:
            rootp = base_dir
        elif self.case_dir:
            rootp = self.case_dir

        model.setRootPath(rootp)

        tree = QtWidgets.QTreeView(None)

        tree.setModel(model)
        tree.setSortingEnabled(True)
        tree.setWindowTitle('Explorer')
        tree.setRootIndex(model.index(rootp))

        nc = tree.header().count()

        for i in range(1, nc):
            tree.hideColumn(i)

        tree.setContextMenuPolicy(QtCore.Qt.ContextMenuPolicy.CustomContextMenu)
        tree.customContextMenuRequested.connect(self.explorerContextMenu)

        tree.doubleClicked.connect(self._explorerDoubleClick)

        return tree;
    # ---------------------------------------------------------------


    # ---------------------------------------------------------------
    def _editSelectedFile(self):
        """
        Edit action for mouse right-click
        """

        self.readOnly = False

        t = "Editor: %s" % (self._currentSelection['filename'])
        self.setWindowTitle(t)

        fn = os.path.join(self.case_dir,
                          self._currentSelection['subpath'],
                          self._currentSelection['filename'])
        self.openFile(fn=fn)
    # ---------------------------------------------------------------


    # ---------------------------------------------------------------
    def _viewSelectedFile(self):
        """
        View action for mouse left-click
        """

        self.readOnly = True

        t = "Viewer: %s" % (self._currentSelection['filename'])
        self.setWindowTitle(t)

        fn = os.path.join(self.case_dir,
                          self._currentSelection['subpath'],
                          self._currentSelection['filename'])
        self.openFile(fn=fn)
    # ---------------------------------------------------------------


    # ---------------------------------------------------------------
    def _copySelectedFile(self):
        """
        Copy files to the SRC folder.
        """

        src_path = os.path.join(self.case_dir,
                                self._currentSelection['subpath'],
                                self._currentSelection['filename'])

        if self.case_name in ('SRC', 'DATA'):
            trg_path = os.path.join(self.case_dir,
                                    self._currentSelection['filename'])
        else:
            sp = self._currentSelection['subpath']
            while '/' in sp and len(sp) > 3:
                e1, e2 = os.path.split(sp)
                if e2 in ('SRC', 'DATA'):
                    break
                else:
                    sp = e1

            trg_path = os.path.join(sp, self._currentSelection['filename'])

        shutil.copy2(src_path, trg_path)
    # ---------------------------------------------------------------


    # ---------------------------------------------------------------
    def _removeSelectedFile(self):
        """
        Remove a file from the SRC dir
        """

        title = "Remove file"
        question = "Remove %s from the SRC folder (Stored in DRAFT) ?" % (self._currentSelection['filename'])

        choice = QMessageBox.question(self, title, question, _msgbox_yes_no())

        if choice == _msgbox_yes():
            fn = os.path.join(self.case_dir,
                              self._currentSelection['subpath'],
                              self._currentSelection['filename'])

            draft = os.path.join(self.case_dir,
                               self._currentSelection['subpath'],
                               'DRAFT')
            if not os.path.exists(draft):
                os.mkdir(draft)
            fn2 = os.path.join(draft, self._currentSelection['filename'])

            if os.path.exists(fn2):
                q = 'A file named %s allready exists in DRAFT.\nDo you want to overwrite it?' % (self._currentSelection['filename'])
                choice2 = QMessageBox.question(self, '', q, _msgbox_yes_no())
                if choice2 == _msgbox_no():
                    return

            shutil.move(fn, fn2)
        else:
            pass
    # ---------------------------------------------------------------


    # ---------------------------------------------------------------
    def _restoreSelectedFile(self):
        """
        Move a file from DRAFT to the SRC folder
        """

        title = "Move to SRC"
        question = "Move file %s from DRAFT to SRC folder ?" % (self._currentSelection['filename'])

        choice = QMessageBox.question(self, title, question, _msgbox_yes_no())

        if choice == _msgbox_yes():
            fn = os.path.join(self.case_dir,
                              self._currentSelection['subpath'],
                              self._currentSelection['filename'])

            fn2 = os.path.join(self.case_dir, self._currentSelection['filename'])

            if os.path.exists(fn2):
                q = 'A file named %s allready exists in SRC\nDo you want to overwrite it?' % (self._currentSelection['filename'])
                choice2 = QMessageBox.question(self, '', q, _msgbox_yes_no())

                if choice2 == _msgbox_no():
                    return

            shutil.move(fn, fn2)

        else:
            pass
    # ---------------------------------------------------------------


    # ---------------------------------------------------------------
    def _deleteSelectedFile(self):
        """
        Remove a file from the SRC dir
        """

        title = "Delete file"
        question = "Really delete %s ?" % (self._currentSelection['filename'])

        choice = QMessageBox.question(self, title, question, _msgbox_yes_no())

        if choice == _msgbox_yes():
            fn = os.path.join(self.case_dir,
                              self._currentSelection['subpath'],
                              self._currentSelection['filename'])

            try:
                os.remove(fn)
            except Exception:
                pass

            d = os.path.split(fn)[0]
            if os.path.basename(d) in ('DRAFT', 'STASH'):
                l = os.listdir(d)
                if len(l) < 1:
                    try:
                        os.rmdir(d)
                    except Exception:
                        pass
        else:
            pass
    # ---------------------------------------------------------------


    # ---------------------------------------------------------------
    def updateFileState(self, new_state = False):
        """
        Update file state (saved or not)
        """
        self.saved  = new_state
        self.textEdit.viewport().update()
    # ---------------------------------------------------------------


    # ---------------------------------------------------------------
    def openFile(self, fn = None):
        """
        Open a file in the editor
        """

        if not self.saved:
            self.closeOpenedFile()

        if fn:
            self.filename = fn
        else:
            self.filename = QFileDialog.getOpenFileName(self, 'Open File',
                                                        self.last_dir)

        if self.filename:
            self.last_dir = os.path.split(self.filename)[0]

        self.textEdit.setReadOnly(self.readOnly)
        self.saveFileAction.setEnabled(not self.readOnly)

        if self.filename != None and self.filename != '':
            file = open(self.filename, 'r')
            self.file_extension = self.filename.split('.')[-1]

            self.newFile()
            with file:
                text = file.read()
                self.textEdit.setPlainText(text)
                self.updateFileState(True)
    # ---------------------------------------------------------------


    # ---------------------------------------------------------------
    def openFileForAction(self, fn = None):

        if self.readOnly != self.readerMode:
            self.readOnly = self.readerMode

        self.openFile(fn)
    # ---------------------------------------------------------------


    # ---------------------------------------------------------------
    def newFile(self):
        """
        Create a new file (blank)
        """

        self.opened = True
        self.updateFileState(False)
        if self.useHighlight:
            hl = QtextHighlighter(self.textEdit.document(), self.file_extension)
        self.textEdit.show()
    # ---------------------------------------------------------------


    # ---------------------------------------------------------------
    def saveFile(self):
        """
        Save file
        """
        if not self.opened:
            return

        if self.filename != None and self.filename != '':
            file = open(self.filename,'w')
            text = self.textEdit.toPlainText()
            file.write(text)
            file.close()

            self.updateFileState(True)

        else:
            self.saveFileAs()
    # ---------------------------------------------------------------


    # ---------------------------------------------------------------
    def saveFileAs(self):
        """
        Save file as
        """
        if not self.opened:
            return

        ret = QFileDialog.getSaveFileName(self, 'Save File')

        if type(ret) == str:
            self.filename = ret
        elif type(ret) == tuple:
            self.filename = ret[0]
        else:
            raise Exception("Uknown return type for 'QFileDialog.getSaveFileName'")

        self.last_dir = os.path.split(self.filename)[0]

        if self.filename != None and self.filename != '':
            file = open(self.filename,'w')
            text = self.textEdit.toPlainText()
            file.write(text)
            file.close()

            self.updateFileState(True)
    # ---------------------------------------------------------------


    # ---------------------------------------------------------------
    def closeOpenedFile(self):
        """
        Close an opened file
        """

        if self.saved == False and self.readOnly == False:
            choice = QMessageBox.question(self, 'Built-in editor',
                                          'File changed.\nDo you want to save?',
                                          _msgbox_yes_no())
            if choice == _msgbox_yes():
                self.saveFile()
            else:
                pass

        self.saved  = True
        self.opened = False

        self.filename = ''
        self.textEdit.setPlainText('')
    # ---------------------------------------------------------------


    # ---------------------------------------------------------------
    def closeApplication(self):
        """
        Close the editor
        """
        if self.opened == True:
            choice = QMessageBox.question(self, 'Built-in editor',
                                          "Exit text editor?",
                                          _msgbox_yes_no())
        else:
            choice = _msgbox_yes()

        if choice == _msgbox_yes():
            self.closeOpenedFile()

            settings = QtCore.QSettings()
            settings.setValue("MainWindow/Geometry",
                              self.saveGeometry())

            self.close()
            return 0
        else:
            return 1
    # ---------------------------------------------------------------


    # ---------------------------------------------------------------
    def closeEvent(self, event):

        decision = self.closeApplication()
        if decision == 1:
            event.ignore()
    # ---------------------------------------------------------------


#-------------------------------------------------------------------------------
# END OF FILE
#-------------------------------------------------------------------------------