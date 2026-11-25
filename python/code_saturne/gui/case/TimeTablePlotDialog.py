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
This module contains the following classes:
- TimeTablePlotDialog
"""

#-------------------------------------------------------------------------------
# Standard modules
#-------------------------------------------------------------------------------

import os
import sys
import logging

#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------

from code_saturne.gui.base.QtCore import *
from code_saturne.gui.base.QtWidgets import *

# Dynamic import supporting PyQt6 and PyQt5 for Matplotlib
try:
    try:
        from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
        from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar
    except ImportError:
        from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
        from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
    from matplotlib.figure import Figure
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

from code_saturne.model.Common import GuiParam

logging.basicConfig()
log = logging.getLogger("TimeTablePlotDialog")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Dialog class
#-------------------------------------------------------------------------------

class TimeTablePlotDialog(QDialog):
    """
    Dialog displaying 2D plots of time table variables as a function of time
    using Matplotlib with variable selection.
    """

    def __init__(self, parent=None, table_name="", times=None, variables_dict=None, time_label="Time (s)"):
        """
        Constructor.

        :param parent: Parent QWidget
        :param table_name: Name of the time table
        :param times: List of time coordinates (float)
        :param variables_dict: Dictionary mapping variable name -> list of float values
        :param time_label: Label for the time axis
        """
        super(TimeTablePlotDialog, self).__init__(parent)

        self.table_name = table_name
        self.times = times if times is not None else []
        self.variables_dict = variables_dict if variables_dict is not None else {}
        self.time_label = time_label

        title = self.tr("Plot Time Table")
        if self.table_name:
            title += f" - {self.table_name}"
        self.setWindowTitle(title)
        self.resize(750, 550)
        self.setMinimumSize(500, 400)

        self._setup_ui()
        self._populate_variables()

        if HAS_MATPLOTLIB and self.variables_dict:
            self.update_plot()


    def _setup_ui(self):
        """
        Set up the dialog layout and widgets.
        """
        main_layout = QVBoxLayout(self)
        main_layout.setContentsMargins(10, 10, 10, 10)
        main_layout.setSpacing(6)

        # Top control bar
        top_bar_layout = QHBoxLayout()
        top_bar_layout.setContentsMargins(0, 0, 0, 0)

        label_var = QLabel(self.tr("Variable:"), self)
        top_bar_layout.addWidget(label_var)

        self.comboBoxVariable = QComboBox(self)
        self.comboBoxVariable.setMinimumWidth(200)
        top_bar_layout.addWidget(self.comboBoxVariable)

        top_bar_layout.addStretch()

        self.buttonClose = QPushButton(self.tr("Close"), self)
        self.buttonClose.clicked.connect(self.accept)
        top_bar_layout.addWidget(self.buttonClose)

        main_layout.addLayout(top_bar_layout)

        # Matplotlib Canvas / Fallback
        if HAS_MATPLOTLIB:
            self.figure = Figure(figsize=(7, 5), dpi=100, tight_layout=True)
            self.canvas = FigureCanvas(self.figure)
            self.canvas.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)

            self.ax = self.figure.add_subplot(111)

            self.toolbar = NavigationToolbar(self.canvas, self)
            main_layout.addWidget(self.toolbar)
            main_layout.addWidget(self.canvas)
        else:
            self.canvas = None
            self.toolbar = None
            self.ax = None

            missing_label = QLabel(
                self.tr("The 'matplotlib' Python package is required to plot time tables,\n"
                        "but it is not installed in the current environment.\n\n"
                        "Please install matplotlib to enable this feature."),
                self
            )
            missing_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
            missing_label.setStyleSheet("color: #cc0000; font-size: 11pt;")
            main_layout.addWidget(missing_label)

        # Bottom status bar
        self.labelStatus = QLabel(self)
        self.labelStatus.setStyleSheet("color: #555555; font-size: 9pt;")
        main_layout.addWidget(self.labelStatus)

        # Connections
        self.comboBoxVariable.currentIndexChanged.connect(self.slotVariableChanged)


    def _populate_variables(self):
        """
        Populate the variable selection combobox.
        """
        self.comboBoxVariable.blockSignals(True)
        self.comboBoxVariable.clear()

        for var_name in self.variables_dict.keys():
            self.comboBoxVariable.addItem(var_name)

        self.comboBoxVariable.blockSignals(False)

        if not self.variables_dict:
            self.comboBoxVariable.setEnabled(False)
            self.labelStatus.setText(self.tr("No variable data available to plot."))


    @Slot(int)
    def slotVariableChanged(self, index):
        """
        Slot called when the selected variable in the combobox changes.
        """
        if HAS_MATPLOTLIB:
            self.update_plot()


    def update_plot(self):
        """
        Update the Matplotlib plot with the currently selected variable.
        """
        if not HAS_MATPLOTLIB or self.ax is None:
            return

        var_name = self.comboBoxVariable.currentText()
        if not var_name or var_name not in self.variables_dict:
            self.ax.clear()
            self.ax.text(
                0.5, 0.5, self.tr("No data to display"),
                horizontalalignment='center', verticalalignment='center',
                transform=self.ax.transAxes
            )
            self.canvas.draw()
            return

        y_values = self.variables_dict[var_name]

        self.ax.clear()

        # Plot line with subtle markers
        self.ax.plot(
            self.times,
            y_values,
            color='#1f77b4',
            linewidth=1.5,
            marker='.',
            markersize=4,
            label=var_name
        )

        self.ax.set_xlabel(self.time_label)
        self.ax.set_ylabel(var_name)
        self.ax.set_title(f"{var_name} vs {self.time_label}")
        self.ax.grid(True, linestyle='--', alpha=0.6)

        self.canvas.draw()

        # Update status message
        n_pts = len(self.times)
        if n_pts > 0 and len(y_values) > 0:
            t_min, t_max = min(self.times), max(self.times)
            y_min, y_max = min(y_values), max(y_values)
            self.labelStatus.setText(
                self.tr("Points: %d  |  Time: [%.4g, %.4g]  |  %s: [%.4g, %.4g]")
                % (n_pts, t_min, t_max, var_name, y_min, y_max)
            )
        else:
            self.labelStatus.setText(self.tr("No points in time table."))
