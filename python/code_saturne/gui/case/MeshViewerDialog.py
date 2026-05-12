# -*- coding: utf-8 -*-

# ============================================================================
# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2026 EDF S.A.
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
# more details.
# ============================================================================

"""
MeshViewerDialog — Mesh visualization window integrated into the GUI
====================================================================

This module provides a 3D Code_Saturne mesh visualization window
based on pyvista/pyvistaqt, integrated into the PyQt GUI.

Optional dependencies (buttons disabled if unavailable):
    pip install pyvista pyvistaqt

Supported mesh format:
    Code_Saturne native .csm file (Big Endian binary format)
    Advantage over .med: directly exposes all boundary faces
    via face_cells[i][1]==0, including polyhedron faces.

Three usage modes
-----------------
1. Boundary face groups mode (SolutionDomainView)
   - Displays colored boundary face groups
   - Click to highlight, automatic legend
   - List of boundary zones if case is provided

2. Boundary zones mode (LocalizationView)
   - OR/AND/NOT expression parser on group names
   - Selected zone highlighted, rest translucent

3. Probes mode (OutputControlView)
   - 3D spheres per probe, red highlight on click
   - Add/remove probes dialog
   - Save to setup.xml + GUI sync via probes_updated signal

Important technical points
--------------------------
- QTimer.singleShot(300, ...): avoids white screen by letting exec_()
  open the window before the first VTK render
- MeshLoaderThread: CSM parsing + subsets/group_pairs/spheres precomputation
  in a background thread to keep GUI responsive
- Pre-extracted subsets and group_pairs: max 2 add_mesh calls on click
- Individual actors per probe: required for add/delete/highlight
- Signal disconnection in closeEvent: prevents VTK crashes on teardown
"""

#  Standard imports
import os
import logging
import numpy as np

#  Qt imports (via Code_Saturne wrappers)
from code_saturne.gui.base.QtCore    import *
from code_saturne.gui.base.QtGui     import *
from code_saturne.gui.base.QtWidgets import *

# QTimer, QThread and pyqtSignal are not always exported by wrappers

if QT_API == "PYSIDE6":
    try:
        import Signal as pyqtSignal # Alias
    except ImportError:
        pass
###try:
###    from PyQt5.QtCore import QTimer, QThread, pyqtSignal
###except ImportError:
###    try:
###        from PySide6.QtCore import QTimer, QThread, Signal as pyqtSignal
###    except ImportError:
###        pass

from code_saturne.model.Common import GuiParam

# ============================================================================
# Logger
# ============================================================================
logging.basicConfig()
log = logging.getLogger("MeshViewerDialog")
log.setLevel(GuiParam.DEBUG)

# ============================================================================
# Imports PyVista modules if available
# ============================================================================
from code_saturne.model.PyVistaWrapper import pyvistaWrapper

_wrapper = pyvistaWrapper()
HAS_VISU = _wrapper.isActive()
if HAS_VISU:
    pv = _wrapper.getPyVista()
    QtInteractor = _wrapper.getQtInteractor()
else:
    pv = None
    QtInteractor = None
    log.warning("MeshViewerDialog: missing pyvista/pyvistaqt dependency, viewer disabled.")

# ============================================================================
# Import ".csm" reader
# ============================================================================
from code_saturne.base.cs_csm_reader import (
    read_csm,
    decode_csm_groups,
    csm_to_pyvista_boundary,
    parse_localization,
    CsmReader,
)

# ============================================================================
# Mesh loading (read) thread
# ============================================================================

class MeshLoaderThread(QThread):
    """
    Asynchronous loader thread for CSM mesh.

    Performs in the background thread (without blocking the GUI):
    1. Binary .csm file parsing
    2. PyVista grid construction
    3. Pre-extraction of subsets and group_pairs per group
    4. Precomputation of pyvista spheres for probes

    Signals
    -------
    finished(grid, gc_to_name, name_to_tags, subsets, group_pairs, probe_spheres)
    error(traceback_str)
    """
    finished = pyqtSignal(object, object, object, object, object, object)
    error    = pyqtSignal(str)

    def __init__(self, csm_file, probes=None):
        """
        Parameters
        ----------
        csm_file : str
            Path to the .csm file
        probes : list of (str, float, float, float) or None
            List of probes (label, x, y, z) for sphere precomputation
        """
        super().__init__()
        self.csm_file = csm_file
        self.probes   = probes

    def run(self):
        try:
            # CSM parsing + subset pre-extraction
            grid, gc_to_name, name_to_tags, subsets, group_pairs = \
                csm_to_pyvista_boundary(self.csm_file)

            # Precompute spheres for probes
            # Radius is proportional to mesh bounding box diagonal (0.5%)
            probe_spheres = []
            if self.probes:
                bounds = grid.bounds
                diag   = ((bounds[1]-bounds[0])**2 +
                          (bounds[3]-bounds[2])**2 +
                          (bounds[5]-bounds[4])**2) ** 0.5
                r = diag * 0.005
                for _, x, y, z in self.probes:
                    probe_spheres.append(pv.Sphere(radius=r, center=(x, y, z)))

            self.finished.emit(grid, gc_to_name, name_to_tags,
                               subsets, group_pairs, probe_spheres)
        except Exception as e:
            import traceback
            self.error.emit(traceback.format_exc())


# ============================================================================
# Add probe dialog
# ============================================================================

class AddProbeDialog(QDialog):
    """
    Input dialog for adding a new monitoring probe.

    Fields: Label (str), X (float), Y (float), Z (float)
    Validation: numeric coordinates, non-empty label
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Ajouter une sonde")
        self.setModal(True)
        self.resize(300, 180)

        layout = QVBoxLayout(self)
        form   = QFormLayout()

        self.lineLabel = QLineEdit()
        self.lineLabel.setPlaceholderText("ex: probe_inlet")
        form.addRow("Label :", self.lineLabel)

        self.lineX = QLineEdit("0.0")
        self.lineY = QLineEdit("0.0")
        self.lineZ = QLineEdit("0.0")
        form.addRow("X :", self.lineX)
        form.addRow("Y :", self.lineY)
        form.addRow("Z :", self.lineZ)

        layout.addLayout(form)

        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(self._on_accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    def _on_accept(self):
        """Validate inputs before accepting."""
        try:
            float(self.lineX.text())
            float(self.lineY.text())
            float(self.lineZ.text())
        except ValueError:
            QMessageBox.warning(self, "Erreur",
                "Les coordonnées X, Y, Z doivent être des nombres.")
            return
        if not self.lineLabel.text().strip():
            QMessageBox.warning(self, "Erreur", "Le label ne peut pas être vide.")
            return
        self.accept()

    def get_probe(self):
        """Return (label, x, y, z) after validation."""
        return (
            self.lineLabel.text().strip(),
            float(self.lineX.text()),
            float(self.lineY.text()),
            float(self.lineZ.text()),
        )


#  Color palette for boundary face groups
COLORS = [
    "#4e9af1",  # blue
    "#f1754e",  # orange
    "#6dbf67",  # green
    "#f1c94e",  # yellow
    "#b36dbf",  # purple
    "#4ecff1",  # cyan
    "#f14e8a",  # pink
    "#9af14e",  # light green
    "#f1a44e",  # light orange
    "#f19a4e",  # salmon
    "#4ef1a4",  # mint
    "#a44ef1",  # mauve
    "#f14e4e",  # red
]


# ============================================================================
# Main visualization dialog
# ============================================================================

class MeshViewerDialog(QDialog):
    """
    3D floating visualization window for Code_Saturne mesh.

    Uses native .csm format for boundary face extraction.
    Three modes based on constructor parameters.

    Signals
    -------
    probes_updated : pyqtSignal
        Emitted after saving probes → triggers refresh in
        OutputControlView._refresh_monitoring_table()

    Constructor Parameters
    ----------------------
    parent : QWidget
        Parent widget (calling GUI view)
    csm_file : str
        Path to the .csm file to visualize
    zone_label : str or None
        Boundary zone label to pre-select (zones mode)
    zone_expr : str or None
        Boundary zone expression (zones mode)
    case : QtCase or None
        Code_Saturne case object (required to load zones and save probes)
    probes : list of (str, float, float, float) or None
        List of probes (label, x, y, z) — activates probes mode
    """

    # Signal emitted after saving probes to XML
    probes_updated = pyqtSignal()

    def __init__(self, parent, csm_file, zone_label=None, zone_expr=None,
                 case=None, probes=None):
        super().__init__(parent)

        #  Attributes
        self.csm_file   = csm_file
        self.zone_label = zone_label
        self.zone_expr  = zone_expr
        self.case       = case
        # Local copy of probes: editable without touching XML until "Apply"
        self.probes     = list(probes) if probes else None

        # Warning flag on close if unsaved modifications exist
        self._probes_modified = False

        #  Window title depending on mode
        if self.probes is not None:
            self.setWindowTitle("Probes — %s" % os.path.basename(csm_file))
        elif zone_label:
            self.setWindowTitle("Zone : %s — %s" % (zone_label,
                                                      os.path.basename(csm_file)))
        else:
            self.setWindowTitle("Mesh viewer — %s" %
                                os.path.basename(csm_file))

        self.resize(1100, 650)
        self.setModal(False)  # non-blocking: main GUI remains interactive

        #  Startup
        if not HAS_VISU:
            self._build_fallback_ui()
            return

        pv.global_theme.allow_empty_mesh = True
        self._build_loading_ui()

        # Start loading in background thread
        self._loader = MeshLoaderThread(csm_file, probes=self.probes)
        self._loader.finished.connect(self._on_mesh_loaded)
        self._loader.error.connect(self._on_mesh_error)
        self._loader.start()

    # =========================================================================
    # Asynchronous loading
    # =========================================================================

    def _build_loading_ui(self):
        """Display a message while loading mesh."""
        layout = QVBoxLayout()
        self.setLayout(layout)
        lbl = QLabel("<b>Currently loading mesh...</b>")
        lbl.setAlignment(Qt.AlignCenter)
        layout.addWidget(lbl)

    def _on_mesh_loaded(self, grid, gc_to_name, name_to_tags,
                        subsets, group_pairs, probe_spheres):
        """
        Callback called when thread has finished loading mesh.
        Stores data and builds UI.
        """
        self.grid           = grid
        self.gc_to_name     = gc_to_name
        self.name_to_tags   = name_to_tags
        self.subsets        = subsets
        self.group_pairs    = group_pairs
        self._probe_spheres = probe_spheres  # spheres precomputed in thread

        # Replace loading layout with main UI
        # (Qt technique to delete an existing layout)
        QWidget().setLayout(self.layout())
        self._build_ui()

    def _on_mesh_error(self, error_msg):
        """Callback on loading error."""
        QWidget().setLayout(self.layout())
        self._build_error_ui(error_msg)

    # =========================================================================
    # UI construction
    # =========================================================================

    def _build_ui(self):
        """
        Build main UI with left panel and 3D viewer.

        Note on QTimer.singleShot:
            PyVista rendering must be triggered AFTER exec_() makes
            the window visible. Without this delay, rendering happens in
            an unmapped widget and the screen remains blank.
        """
        layout = QHBoxLayout(self)
        left   = QVBoxLayout()

        # Left panel depending on mode
        if self.probes is not None:
            # Probes mode: initialize _unique_names for mesh rendering
            self._unique_names = list(self.name_to_tags.keys())
            self._build_probe_panel(left)
        else:
            # Groups/zones mode
            self._build_group_panel(left)

        # Color map indexed by group name
        self._color_map = {name: COLORS[i % len(COLORS)]
                           for i, name in enumerate(self._unique_names)}

        left_widget = QWidget()
        left_widget.setLayout(left)
        left_widget.setFixedWidth(270)
        layout.addWidget(left_widget)

        # VTK viewer (QtInteractor = standard QWidget)
        self.plotter = QtInteractor(self)
        layout.addWidget(self.plotter.interactor)

        # Viewer internal state
        self._first_render       = True   # True → view_isometric() on first render
        self._probe_actors       = {}     # {index: pyvista sphere actor}
        self._selected_probe_idx = None   # selected probe index

        # Trigger first render after window is displayed
        if self.probes is not None:
            QTimer.singleShot(300, self._show_probes)
        elif self.zone_expr:
            QTimer.singleShot(300,
                lambda: self._show_zone(self.zone_label, self.zone_expr))
        else:
            QTimer.singleShot(300, self._show_all)

    def _build_group_panel(self, left):
        """Build left panel in boundary face groups / boundary zones mode."""
        left.addWidget(QLabel("<b>Boundary faces' groups</b>"))

        self.group_list = QListWidget()
        self._unique_names = list(self.name_to_tags.keys())
        for name in self._unique_names:
            self.group_list.addItem(name)
        self.group_list.currentRowChanged.connect(self._on_group_selected)
        left.addWidget(self.group_list)

        btn_reset = QPushButton("Show all")
        btn_reset.clicked.connect(self._show_all)
        left.addWidget(btn_reset)

        btn_cam = QPushButton("Reset camera")
        btn_cam.clicked.connect(self._reset_camera)
        left.addWidget(btn_cam)

        # Boundary zones list (only if case provided)
        if self.case is not None:
            separator = QFrame()
            separator.setFrameShape(QFrame.HLine)
            separator.setFrameShadow(QFrame.Sunken)
            left.addWidget(separator)

            left.addWidget(QLabel("<b>Boundary Zones</b>"))

            self.zone_list = QListWidget()
            self._load_zones()
            self.zone_list.currentRowChanged.connect(self._on_zone_selected)
            left.addWidget(self.zone_list)

    def _build_probe_panel(self, left):
        """Build left panel in probes mode."""
        left.addWidget(QLabel("<b>Probes</b>"))

        self.probe_list = QListWidget()
        self._refresh_probe_list()
        self.probe_list.currentRowChanged.connect(self._on_probe_selected)
        left.addWidget(self.probe_list)

        # Actions on probes
        btn_add = QPushButton("Add a probe")
        btn_add.clicked.connect(self._on_add_probe)
        left.addWidget(btn_add)

        btn_del = QPushButton("Delete probe")
        btn_del.clicked.connect(self._on_delete_probe)
        left.addWidget(btn_del)

        # Visual separator
        separator = QFrame()
        separator.setFrameShape(QFrame.HLine)
        separator.setFrameShadow(QFrame.Sunken)
        left.addWidget(separator)

        # Apply button (save to XML)
        btn_apply = QPushButton("Apply")
        btn_apply.setToolTip("Save modifications to setup.xml")
        btn_apply.setStyleSheet(
            "QPushButton { background-color: #4e9af1; color: white; "
            "font-weight: bold; padding: 4px; }")
        btn_apply.clicked.connect(self._on_apply_probes)
        left.addWidget(btn_apply)

        btn_reset = QPushButton("Show all")
        btn_reset.clicked.connect(self._show_probes)
        left.addWidget(btn_reset)

        btn_cam = QPushButton("Reset camera")
        btn_cam.clicked.connect(self._reset_camera)
        left.addWidget(btn_cam)

    def _build_fallback_ui(self):
        """Fallback UI displayed if pyvista/pyvistaqt are not installed."""
        layout = QVBoxLayout()
        self.setLayout(layout)
        msg = QLabel(
            "<b>Viewer not available.</b><br><br>"
            "The following libraries are required:<br>"
            "<code>pip install pyvista pyvistaqt</code><br><br>"
        )
        msg.setWordWrap(True)
        layout.addWidget(msg)

    def _build_error_ui(self, error_msg):
        """UI displayed on mesh loading error."""
        layout = QVBoxLayout()
        self.setLayout(layout)
        layout.addWidget(QLabel("<b>Error during mesh loading step :</b>"))
        txt = QTextEdit()
        txt.setPlainText(error_msg)
        txt.setReadOnly(True)
        layout.addWidget(txt)

    def _load_zones(self):
        """Load boundary zones from setup.xml via LocalizationModel."""
        try:
            from code_saturne.model.LocalizationModel import LocalizationModel
            mdl = LocalizationModel("BoundaryZone", self.case)
            self._zones = [(z.getLabel(), z.getLocalization())
                           for z in mdl.getZones()]
            for label, expr in self._zones:
                self.zone_list.addItem("%s  —  %s" % (label, expr))
        except Exception as e:
            log.warning("_load_zones error: %s" % e)
            self._zones = []

    # =========================================================================
    # Camera management
    # =========================================================================

    def _get_camera_state(self):
        """
        Save current camera state.
        Returns None if plotter is not yet initialized.
        """
        try:
            return {
                'position':    self.plotter.camera.position,
                'focal_point': self.plotter.camera.focal_point,
                'view_up':     self.plotter.camera.up,
            }
        except Exception:
            return None

    def _set_camera_state(self, state):
        """Restore camera state saved by _get_camera_state."""
        if state is None:
            return
        try:
            self.plotter.camera.position    = state['position']
            self.plotter.camera.focal_point = state['focal_point']
            self.plotter.camera.up          = state['view_up']
        except Exception:
            pass

    def _reset_camera(self):
        """Reset camera to isometric view (connected to reset button)."""
        if not hasattr(self, 'plotter'):
            return
        self.plotter.view_isometric()
        self.plotter.render()

    def _probe_radius(self):
        """
        Compute probe sphere radius.
        Proportional to mesh bounding box diagonal (0.5%) to automatically
        adapt to geometry scale.
        """
        try:
            bounds = self.grid.bounds
            diag = ((bounds[1]-bounds[0])**2 +
                    (bounds[3]-bounds[2])**2 +
                    (bounds[5]-bounds[4])**2) ** 0.5
            return diag * 0.005
        except Exception:
            return 0.01

    def _finalize_render(self, cam):
        """
        Finalize render: reset clipping range, camera, then call render().

        On first render (_first_render=True), forces view_isometric().
        Afterwards, restores saved camera state to preserve viewpoint.
        """
        self.plotter.renderer.ResetCameraClippingRange()
        if self._first_render or not cam:
            self.plotter.view_isometric()
            self._first_render = False
        else:
            self._set_camera_state(cam)
        self.plotter.render()

    # =========================================================================
    # Boundary face groups mode
    # =========================================================================

    def _show_all(self):
        """
        Display all face groups in solid colors with legend.
        Uses subsets pre-extracted in thread (1 add_mesh call per group).
        """
        if not hasattr(self, 'plotter'):
            return

        if hasattr(self, 'group_list'):
            self.group_list.clearSelection()

        cam = self._get_camera_state()
        self.plotter.clear()

        for name, subset in self.subsets.items():
            self.plotter.add_mesh(subset,
                                   color=self._color_map[name],
                                   opacity=1.0,
                                   show_edges=True,
                                   edge_color="gray",
                                   reset_camera=False)

        legend_entries = [[n, self._color_map[n]] for n in self.subsets]
        if legend_entries:
            self.plotter.add_legend(legend_entries)

        self._finalize_render(cam)

    def _on_group_selected(self, row):
        """
        Highlight selected group in the list.

        Uses group_pairs pre-extracted in thread:
        - group faces: solid color, black edges
        - remaining faces: translucent gray (opacity 0.3)
        → Only 2 add_mesh calls, no extract_cells on click.
        """
        if not hasattr(self, 'plotter'):
            return
        if row < 0:
            return
        if hasattr(self, 'zone_list'):
            self.zone_list.clearSelection()

        cam  = self._get_camera_state()
        name = self._unique_names[row]
        self.plotter.clear()

        if name in self.group_pairs:
            sel_grid, rest_grid = self.group_pairs[name]
            self.plotter.add_mesh(sel_grid,
                                   color=self._color_map[name],
                                   opacity=1.0,
                                   show_edges=True,
                                   edge_color="black",
                                   reset_camera=False)
            if rest_grid.n_cells > 0:
                self.plotter.add_mesh(rest_grid,
                                       color="lightgray",
                                       opacity=0.3,
                                       show_edges=False,
                                       reset_camera=False)

        self.plotter.add_legend([[name, self._color_map[name]]])
        self._finalize_render(cam)

    # =========================================================================
    # Boundary zones mode
    # =========================================================================

    def _on_zone_selected(self, row):
        """Callback when user selects a boundary zone."""
        if not hasattr(self, 'plotter'):
            return
        if row < 0:
            return
        self.group_list.clearSelection()
        label, expr = self._zones[row]
        self._show_zone(label, expr)

    def _show_zone(self, label, expr):
        """
        Display a boundary zone identified by its expression.

        Parses the expression (OR/AND/NOT) and displays:
        - zone faces: blue (#4e9af1), black edges
        - remaining faces: translucent gray (opacity 0.3)
        """
        if not hasattr(self, 'plotter'):
            return

        try:
            cam = self._get_camera_state()
            self.plotter.clear()

            mask = parse_localization(expr, self.grid, self.name_to_tags)

            if mask is None or not mask.any():
                self.plotter.add_text(
                    "No face found for : %s" % expr, font_size=12)
                self.plotter.render()
                return

            zone_grid = self.grid.extract_cells(mask)
            rest_grid = self.grid.extract_cells(~mask)

            if zone_grid.n_cells > 0:
                self.plotter.add_mesh(zone_grid,
                                       color="#4e9af1",
                                       opacity=1.0,
                                       show_edges=True,
                                       edge_color="black",
                                       reset_camera=False)
            if rest_grid.n_cells > 0:
                self.plotter.add_mesh(rest_grid,
                                       color="lightgray",
                                       opacity=0.3,
                                       show_edges=False,
                                       reset_camera=False)
            if zone_grid.n_cells > 0:
                self.plotter.add_legend([[label, "#4e9af1"]])

            self._finalize_render(cam)

        except Exception as e:
            import traceback
            traceback.print_exc()

    # =========================================================================
    # Probes mode
    # =========================================================================

    def _refresh_probe_list(self):
        """Update QListWidget from self.probes (local list)."""
        self.probe_list.clear()
        for name, x, y, z in self.probes:
            self.probe_list.addItem(
                "%s  (%.3f, %.3f, %.3f)" % (name, x, y, z))

    def _show_probes(self):
        """
        Display translucent mesh + all probes as spheres.

        Uses spheres precomputed in thread (_probe_spheres).
        Creates an individual actor per probe to support add/delete/highlight.
        """
        if not hasattr(self, 'plotter'):
            return

        cam = self._get_camera_state()
        self.plotter.clear()
        self._probe_actors       = {}
        self._selected_probe_idx = None

        # Translucent mesh (geometric background)
        for name, subset in self.subsets.items():
            self.plotter.add_mesh(subset,
                                   color="lightgray",
                                   opacity=0.2,
                                   show_edges=False,
                                   reset_camera=False)

        # Individual spheres per probe
        r = self._probe_radius()
        for i, (name, x, y, z) in enumerate(self.probes):
            # Use precomputed sphere if available
            if i < len(self._probe_spheres):
                sphere = self._probe_spheres[i]
            else:
                sphere = pv.Sphere(radius=r, center=(x, y, z))
            actor = self.plotter.add_mesh(sphere, color="#4e9af1",
                                           reset_camera=False)
            self._probe_actors[i] = actor

        self._finalize_render(cam)

    def _on_probe_selected(self, row):
        """
        Highlight selected probe (red sphere 2x larger).

        Since existing VTK actor geometry cannot be modified in place,
        remove and recreate the relevant spheres (previous selection → small blue,
        new selection → large red).
        """
        if not hasattr(self, 'plotter'):
            return
        if row < 0:
            return
        if not self._probe_actors:
            self._show_probes()
            return

        cam = self._get_camera_state()
        r   = self._probe_radius()

        # Restore previous selected probe as small blue sphere
        if self._selected_probe_idx is not None:
            old_idx = self._selected_probe_idx
            if old_idx in self._probe_actors:
                name, x, y, z = self.probes[old_idx]
                self.plotter.remove_actor(self._probe_actors[old_idx])
                sphere = (self._probe_spheres[old_idx]
                          if old_idx < len(self._probe_spheres)
                          else pv.Sphere(radius=r, center=(x, y, z)))
                actor = self.plotter.add_mesh(sphere, color="#4e9af1",
                                               reset_camera=False)
                self._probe_actors[old_idx] = actor

        # Replace newly selected probe with large red sphere
        if row < len(self.probes):
            name, x, y, z = self.probes[row]
            self.plotter.remove_actor(self._probe_actors[row])
            big_sphere = pv.Sphere(radius=r * 2.0, center=(x, y, z))
            actor = self.plotter.add_mesh(big_sphere, color="red",
                                           reset_camera=False)
            self._probe_actors[row]  = actor
            self._selected_probe_idx = row

        self.plotter.renderer.ResetCameraClippingRange()
        if cam:
            self._set_camera_state(cam)
        self.plotter.render()

    def _on_add_probe(self):
        """
        Open AddProbeDialog and add probe to local list.
        Updates plotter in real time (no full redraw).
        """
        dlg = AddProbeDialog(self)
        if dlg.exec_() != QDialog.Accepted:
            return

        label, x, y, z = dlg.get_probe()
        self.probes.append((label, x, y, z))

        # Add sphere to plotter
        r      = self._probe_radius()
        sphere = pv.Sphere(radius=r, center=(x, y, z))
        self._probe_spheres.append(sphere)
        actor = self.plotter.add_mesh(sphere, color="#4e9af1",
                                       reset_camera=False)
        self._probe_actors[len(self.probes) - 1] = actor

        self._refresh_probe_list()
        self._probes_modified = True
        self.plotter.renderer.ResetCameraClippingRange()
        self.plotter.render()

    def _on_delete_probe(self):
        """
        Delete selected probe after confirmation.

        Rebuilds _probe_actors dictionary so indices remain consecutive
        without gaps after deletion.
        """
        row = self.probe_list.currentRow()
        if row < 0:
            QMessageBox.information(self, "Delete",
                "Please select a probe from the list.")
            return

        reply = QMessageBox.question(self,
            "Confirm probe deletion",
            "Delete probe '%s' ?" % self.probes[row][0],
            QMessageBox.Yes | QMessageBox.No)
        if reply != QMessageBox.Yes:
            return

        # Delete VTK actor
        if row in self._probe_actors:
            self.plotter.remove_actor(self._probe_actors[row])

        # Rebuild dictionary shifting indices > row by -1
        new_actors = {}
        for i, actor in self._probe_actors.items():
            if i < row:
                new_actors[i] = actor
            elif i > row:
                new_actors[i - 1] = actor
        self._probe_actors = new_actors

        # Delete corresponding precomputed sphere
        if row < len(self._probe_spheres):
            self._probe_spheres.pop(row)

        # Remove probe from local list
        self.probes.pop(row)

        self._selected_probe_idx = None
        self._refresh_probe_list()
        self._probes_modified = True
        self.plotter.renderer.ResetCameraClippingRange()
        self.plotter.render()

    def _on_apply_probes(self):
        """
        Save probes from self.probes to setup.xml.

        Strategy: delete all existing probes then recreate from
        local list (ensures consecutive numbering without gaps).

        Emits probes_updated to trigger refresh in OutputControlView.
        """
        if self.case is None:
            QMessageBox.warning(self, "Error",
                "Impossible to save : case not available.")
            return

        try:
            from code_saturne.model.OutputControlModel import OutputControlModel
            mdl = OutputControlModel(self.case)

            # Delete all existing probes (from last to first)
            n = mdl.getNumberOfMonitoringPoints()
            for i in range(n, 0, -1):
                mdl.deleteMonitoringPoint(str(i))

            # Recreate probes from local list with 1..N numbering
            for i, (label, x, y, z) in enumerate(self.probes):
                mdl.addMonitoringPoint(x, y, z)
                num = str(i + 1)
                try:
                    mdl.setMonitoringPointName(num, label)
                except Exception:
                    pass  # setMonitoringPointName might not exist

            self._probes_modified = False
            # Notify OutputControlView to refresh its table
            self.probes_updated.emit()
            QMessageBox.information(self, "Saved",
                "%d probes saved to setup.xml." % len(self.probes))

        except Exception as e:
            import traceback
            traceback.print_exc()
            QMessageBox.critical(self, "Erreur de sauvegarde", str(e))

    # =========================================================================
    # Teardown / Close
    # =========================================================================

    def closeEvent(self, event):
        """
        Handle clean dialog closure.

        1. Warn if probe modifications are unsaved
        2. Disconnect Qt signals to prevent calls to plotter.clear()
           after VTK teardown has started (→ AttributeError)
        3. Wait for loading thread if still running
        4. Cleanly close VTK plotter
        """
        # Warn if unsaved modifications
        if self._probes_modified:
            reply = QMessageBox.question(
                self,
                "Modifications not saved",
                "Some probes related modifications were not saved.\n"
                "Quit without saving ?",
                QMessageBox.Yes | QMessageBox.No)
            if reply == QMessageBox.No:
                event.ignore()
                return

        # Disconnect all list signals to prevent currentRowChanged
        # from being triggered during dialog destruction and calling
        # plotter.clear() on a partially destroyed renderer
        try:
            if hasattr(self, 'group_list'):
                self.group_list.currentRowChanged.disconnect()
            if hasattr(self, 'zone_list'):
                self.zone_list.currentRowChanged.disconnect()
            if hasattr(self, 'probe_list'):
                self.probe_list.currentRowChanged.disconnect()
        except Exception:
            pass  # disconnect might fail if signal was not connected

        # Wait for loading thread if still active
        if hasattr(self, '_loader') and self._loader.isRunning():
            self._loader.wait()

        # Cleanly close VTK plotter
        if hasattr(self, 'plotter'):
            self.plotter.close()

        event.accept()
