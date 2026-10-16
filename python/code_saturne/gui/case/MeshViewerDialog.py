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
import math
import re
import logging
import unittest
import numpy as np

#  Qt imports (via Code_Saturne wrappers)
from code_saturne.gui.base.QtCore    import *
from code_saturne.gui.base.QtGui     import *
from code_saturne.gui.base.QtWidgets import *

if QT_API == "PYSIDE6":
    try:
        import Signal as pyqtSignal # Alias
    except ImportError:
        pass

from code_saturne.model.Common import GuiParam

# ============================================================================
# Profile formula evaluation
# ============================================================================

def eval_profile_formula(formula, nb_points, notebook_vars=None):
    """
    Evaluate mathematical formula for 1D profile where parameter s in [0, 1].

    Parameters
    ----------
    formula : str
        Formula defining x, y, z as a function of s.
        Supports C/Python syntax, semicolons, comments (# or //), and power (^).
    nb_points : int, optional
        Number of points along the segment, default 200 (minimum 2).
    notebook_vars : dict, optional
        Dictionary of notebook variables/constants.

    Returns
    -------
    np.ndarray of shape (nb_points, 3)
        Array of (x, y, z) coordinates.
    """
    if not formula or not str(formula).strip():
        raise ValueError("Formula is empty.")

    # 1. Clean comments (// style)
    cleaned = re.sub(r'//.*', '', str(formula))
    # 2. Replace ^ with ** (power)
    cleaned = re.sub(r'\^', '**', cleaned)
    # 3. Replace semicolons with newlines
    cleaned = cleaned.replace(';', '\n')
    # 4. Strip indentation from each line to avoid unexpected indent
    lines = [line.strip() for line in cleaned.splitlines() if line.strip()]
    cleaned = '\n'.join(lines)

    # Base math environment
    math_env = {k: getattr(math, k) for k in dir(math) if not k.startswith('_')}
    math_env['pi'] = math.pi
    math_env['e'] = math.e
    math_env['abs'] = abs
    math_env['min'] = min
    math_env['max'] = max

    if notebook_vars:
        for k, v in notebook_vars.items():
            try:
                math_env[str(k)] = float(v)
            except (ValueError, TypeError):
                pass

    nb_points = max(int(nb_points), 2)
    s_vals = np.linspace(0.0, 1.0, nb_points)

    compiled = compile(cleaned, '<profile_formula>', 'exec')
    points = []
    for s in s_vals:
        scope = dict(math_env)
        scope['s'] = float(s)
        exec(compiled, {'__builtins__': {}}, scope)
        if 'x' not in scope or 'y' not in scope or 'z' not in scope:
            missing = [var for var in ['x', 'y', 'z'] if var not in scope]
            raise ValueError(f"Formula must define variables x, y, and z. Missing: {missing}")
        points.append([float(scope['x']), float(scope['y']), float(scope['z'])])

    return np.array(points, dtype=np.float64)

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

    def __init__(self, parent, csm_file=None, zone_label=None, zone_expr=None,
                 case=None, probes=None, profile=None, profiles=None,
                 profile_display_mode="line", display_mode=None,
                 show_wireframe=True, show_volume=True,
                 stl_file=None, stl_label=None, seed_points=None):

        super().__init__(parent)

        #  Attributes
        self.csm_file   = csm_file
        self.zone_label = zone_label
        self.zone_expr  = zone_expr
        self.case       = case
        # Local copy of probes: editable without touching XML until "Apply"
        self.probes     = list(probes) if probes else None

        # Local copy of profiles (for 1D profile visualization mode)
        if profiles:
            self.profiles = list(profiles)
        elif profile:
            self.profiles = [profile]
        else:
            self.profiles = None

        self._selected_profile_idx = None
        # Profile display mode: 'line' (straight line / tube) or 'points' (discrete probe spheres)
        self.profile_display_mode = display_mode or profile_display_mode or "line"

        # STL mode attributes
        self.stl_file = stl_file
        self.stl_label = stl_label or (os.path.basename(stl_file) if stl_file else None)
        self.seed_points = list(seed_points) if seed_points else []
        self.stl_polydata = None
        self.show_stl = True
        self.show_stl_edges = False
        self.show_seed_points = True
        self.stl_representation = "surface"

        # Wireframe & translucent volume display flags for background mesh
        self.show_wireframe = bool(show_wireframe)
        self.show_volume    = bool(show_volume)

        # Warning flag on close if unsaved modifications exist
        self._probes_modified = False

        #  Window title depending on mode
        mesh_name = os.path.basename(csm_file) if csm_file else "No mesh"
        if self.stl_file is not None:
            lbl = self.stl_label if self.stl_label else "STL"
            self.setWindowTitle("STL Object : %s — %s" % (lbl, mesh_name))
        elif self.probes is not None:
            self.setWindowTitle("Probes — %s" % mesh_name)
        elif self.profiles is not None:
            lbl = self.profiles[0].get('label', '') if len(self.profiles) == 1 else ''
            if lbl:
                self.setWindowTitle("Profile : %s — %s" % (lbl, mesh_name))
            else:
                self.setWindowTitle("Profiles — %s" % mesh_name)
        elif zone_label:
            self.setWindowTitle("Zone : %s — %s" % (zone_label, mesh_name))
        else:
            self.setWindowTitle("Mesh viewer — %s" % mesh_name)

        self.resize(1100, 650)
        self.setModal(False)  # non-blocking: main GUI remains interactive

        #  Startup
        if not HAS_VISU:
            self._build_fallback_ui()
            return

        pv.global_theme.allow_empty_mesh = True

        if csm_file:
            self._build_loading_ui()
            # Start loading in background thread
            self._loader = MeshLoaderThread(csm_file, probes=self.probes)
            self._loader.finished.connect(self._on_mesh_loaded)
            self._loader.error.connect(self._on_mesh_error)
            self._loader.start()
        else:
            self.grid = None
            self.gc_to_name = {}
            self.name_to_tags = {}
            self.subsets = {}
            self.group_pairs = {}
            self._probe_spheres = []
            self._build_ui()

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

        if not hasattr(self, 'show_wireframe'):
            self.show_wireframe = True
        if not hasattr(self, 'show_volume'):
            self.show_volume = True

        # Left panel depending on mode
        if getattr(self, 'stl_file', None) is not None:
            self._unique_names = list(self.name_to_tags.keys()) if getattr(self, 'name_to_tags', None) else []
            self._build_stl_panel(left)
        elif getattr(self, 'probes', None) is not None:
            # Probes mode: initialize _unique_names for mesh rendering
            self._unique_names = list(self.name_to_tags.keys())
            self._build_probe_panel(left)
        elif getattr(self, 'profiles', None) is not None:
            # Profiles mode: initialize _unique_names for mesh rendering
            self._unique_names = list(self.name_to_tags.keys())
            self._build_profile_panel(left)
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
        self.plotter.show_axes()

        # Viewer internal state
        self._first_render         = True   # True → view_isometric() on first render
        self._probe_actors         = {}     # {index: pyvista sphere actor}
        self._selected_probe_idx   = None   # selected probe index
        self._profile_actors       = {}     # {index: (tube, start_sph, end_sph)}
        self._selected_profile_idx = None   # selected profile index

        # Trigger first render after window is displayed
        if getattr(self, 'stl_file', None) is not None:
            QTimer.singleShot(300, self._show_stl)
        elif getattr(self, 'probes', None) is not None:
            QTimer.singleShot(300, self._show_probes)
        elif getattr(self, 'profiles', None) is not None:
            QTimer.singleShot(300, self._show_profiles)
        elif getattr(self, 'zone_expr', None):
            QTimer.singleShot(300,
                lambda: self._show_zone(self.zone_label, self.zone_expr))
        else:
            QTimer.singleShot(300, self._show_all)

    def _build_stl_panel(self, left):
        """Build left panel in STL visualization mode."""
        left.addWidget(QLabel("<b>STL Object Details</b>"))

        if not hasattr(self, 'stl_polydata'):
            self.stl_polydata = None
        if not hasattr(self, 'stl_representation'):
            self.stl_representation = "surface"
        if not hasattr(self, 'show_stl'):
            self.show_stl = True
        if not hasattr(self, 'show_stl_edges'):
            self.show_stl_edges = False
        if not hasattr(self, 'show_seed_points'):
            self.show_seed_points = True
        if not hasattr(self, 'seed_points'):
            self.seed_points = []

        # Load STL geometry info if not already loaded
        stl_file = getattr(self, 'stl_file', None)
        if self.stl_polydata is None and stl_file and os.path.isfile(stl_file):
            try:
                self.stl_polydata = pv.read(stl_file)
            except Exception as e:
                log.warning("Failed to load STL file %s: %s", stl_file, e)

        # Details display
        self.stl_details = QTextEdit()
        self.stl_details.setReadOnly(True)
        self.stl_details.setMaximumHeight(160)

        stl_info = []
        lbl = getattr(self, 'stl_label', None)
        if lbl:
            stl_info.append("<b>Label:</b> %s" % lbl)
        if stl_file:
            stl_info.append("<b>File:</b> %s" % os.path.basename(stl_file))

        stl_bounds = None
        if self.stl_polydata is not None:
            stl_bounds = self.stl_polydata.bounds
            stl_info.append("<b>Points:</b> %d" % self.stl_polydata.n_points)
            stl_info.append("<b>Triangles:</b> %d" % self.stl_polydata.n_cells)
            dx = stl_bounds[1] - stl_bounds[0]
            dy = stl_bounds[3] - stl_bounds[2]
            dz = stl_bounds[5] - stl_bounds[4]
            stl_info.append("<b>Dimensions:</b><br/>&nbsp;ΔX=%.4g, ΔY=%.4g, ΔZ=%.4g" % (dx, dy, dz))
            stl_info.append("<b>Bounds:</b><br/>"
                           "&nbsp;X: [%.4g, %.4g]<br/>"
                           "&nbsp;Y: [%.4g, %.4g]<br/>"
                           "&nbsp;Z: [%.4g, %.4g]" % (
                               stl_bounds[0], stl_bounds[1],
                               stl_bounds[2], stl_bounds[3],
                               stl_bounds[4], stl_bounds[5]))

        self.stl_details.setHtml("<br/>".join(stl_info))
        left.addWidget(self.stl_details)

        # Position check against CFD mesh bounds
        mesh_bounds = None
        if getattr(self, 'grid', None) is not None and hasattr(self.grid, 'bounds'):
            mesh_bounds = self.grid.bounds
        elif getattr(self, 'subsets', None):
            all_pts = []
            for sub in self.subsets.values():
                if hasattr(sub, 'points') and len(sub.points) > 0:
                    all_pts.append(sub.points)
            if all_pts:
                pts_arr = np.vstack(all_pts)
                mesh_bounds = [pts_arr[:, 0].min(), pts_arr[:, 0].max(),
                               pts_arr[:, 1].min(), pts_arr[:, 1].max(),
                               pts_arr[:, 2].min(), pts_arr[:, 2].max()]

        if stl_bounds is not None and mesh_bounds is not None:
            tol = 1e-6
            is_inside = (mesh_bounds[0] - tol <= stl_bounds[0] and
                         stl_bounds[1] <= mesh_bounds[1] + tol and
                         mesh_bounds[2] - tol <= stl_bounds[2] and
                         stl_bounds[3] <= mesh_bounds[3] + tol and
                         mesh_bounds[4] - tol <= stl_bounds[4] and
                         stl_bounds[5] <= mesh_bounds[5] + tol)
            if is_inside:
                self.label_position_check = QLabel("<b style='color:#27ae60;'>✔ STL is INSIDE mesh bounds</b>")
            else:
                self.label_position_check = QLabel("<b style='color:#e67e22;'>⚠ STL extends OUTSIDE mesh bounds</b>")
            self.label_position_check.setAlignment(Qt.AlignCenter)
            left.addWidget(self.label_position_check)

        # Visual separator
        sep = QFrame()
        sep.setFrameShape(QFrame.HLine)
        sep.setFrameShadow(QFrame.Sunken)
        left.addWidget(sep)

        # STL display controls
        left.addWidget(QLabel("<b>STL Display Options</b>"))
        self.check_stl = QCheckBox("Show STL object")
        self.check_stl.setChecked(self.show_stl)
        self.check_stl.toggled.connect(self._on_stl_toggled)
        left.addWidget(self.check_stl)

        left.addWidget(QLabel("Representation:"))
        self.combo_stl_rep = QComboBox()
        self.combo_stl_rep.addItem("Surface", "surface")
        self.combo_stl_rep.addItem("Surface with edges", "surface_edges")
        self.combo_stl_rep.addItem("Wireframe", "wireframe")
        self.combo_stl_rep.addItem("Points", "points")
        self.combo_stl_rep.currentIndexChanged.connect(self._on_stl_rep_changed)
        left.addWidget(self.combo_stl_rep)

        # Visual separator
        sep2 = QFrame()
        sep2.setFrameShape(QFrame.HLine)
        sep2.setFrameShadow(QFrame.Sunken)
        left.addWidget(sep2)

        # Background mesh options
        left.addWidget(QLabel("<b>Background Mesh Options</b>"))
        self.check_wireframe = QCheckBox("Show wireframe")
        self.check_wireframe.setChecked(self.show_wireframe)
        self.check_wireframe.toggled.connect(self._on_wireframe_toggled)
        left.addWidget(self.check_wireframe)

        self.check_volume = QCheckBox("Translucent volume")
        self.check_volume.setChecked(self.show_volume)
        self.check_volume.toggled.connect(self._on_volume_toggled)
        left.addWidget(self.check_volume)

        if self.seed_points:
            sep3 = QFrame()
            sep3.setFrameShape(QFrame.HLine)
            sep3.setFrameShadow(QFrame.Sunken)
            left.addWidget(sep3)

            left.addWidget(QLabel("<b>Seed Points</b>"))
            self.check_seeds = QCheckBox("Show seed points (%d)" % len(self.seed_points))
            self.check_seeds.setChecked(self.show_seed_points)
            self.check_seeds.toggled.connect(self._on_seeds_toggled)
            left.addWidget(self.check_seeds)

        # Button reset camera
        btn_cam = QPushButton("Reset camera")
        btn_cam.clicked.connect(self._reset_camera)
        left.addWidget(btn_cam)

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

        # Visual separator
        separator2 = QFrame()
        separator2.setFrameShape(QFrame.HLine)
        separator2.setFrameShadow(QFrame.Sunken)
        left.addWidget(separator2)

        # Display options
        left.addWidget(QLabel("<b>Display options</b>"))
        self.check_wireframe = QCheckBox("Show wireframe")
        self.check_wireframe.setChecked(self.show_wireframe)
        self.check_wireframe.toggled.connect(self._on_wireframe_toggled)
        left.addWidget(self.check_wireframe)

        self.check_volume = QCheckBox("Translucent volume")
        self.check_volume.setChecked(self.show_volume)
        self.check_volume.toggled.connect(self._on_volume_toggled)
        left.addWidget(self.check_volume)

        btn_reset = QPushButton("Show all")
        btn_reset.clicked.connect(self._show_probes)
        left.addWidget(btn_reset)

        btn_cam = QPushButton("Reset camera")
        btn_cam.clicked.connect(self._reset_camera)
        left.addWidget(btn_cam)

    def _build_profile_panel(self, left):
        """Build left panel in profiles mode."""
        left.addWidget(QLabel("<b>Profiles</b>"))

        self.profile_list = QListWidget()
        for prof in self.profiles:
            lbl = prof.get('label', 'profile')
            self.profile_list.addItem(lbl)
        self.profile_list.currentRowChanged.connect(self._on_profile_selected)
        left.addWidget(self.profile_list)

        # Profile details panel
        left.addWidget(QLabel("<b>Profile details</b>"))
        self.profile_details = QTextEdit()
        self.profile_details.setReadOnly(True)
        self.profile_details.setMaximumHeight(160)
        left.addWidget(self.profile_details)

        # Profile display mode selector
        left.addWidget(QLabel("<b>Display mode</b>"))
        self.combo_profile_display = QComboBox()
        self.combo_profile_display.addItem("Straight line", "line")
        self.combo_profile_display.addItem("Points", "points")
        cur_mode = getattr(self, 'profile_display_mode', 'line')
        idx = 1 if cur_mode in ("points", "spheres") else 0
        self.combo_profile_display.setCurrentIndex(idx)
        self.combo_profile_display.currentIndexChanged.connect(self._on_profile_display_mode_changed)
        left.addWidget(self.combo_profile_display)

        # Visual separator
        separator2 = QFrame()
        separator2.setFrameShape(QFrame.HLine)
        separator2.setFrameShadow(QFrame.Sunken)
        left.addWidget(separator2)

        # Display options
        left.addWidget(QLabel("<b>Display options</b>"))
        self.check_wireframe = QCheckBox("Show wireframe")
        self.check_wireframe.setChecked(self.show_wireframe)
        self.check_wireframe.toggled.connect(self._on_wireframe_toggled)
        left.addWidget(self.check_wireframe)

        self.check_volume = QCheckBox("Translucent volume")
        self.check_volume.setChecked(self.show_volume)
        self.check_volume.toggled.connect(self._on_volume_toggled)
        left.addWidget(self.check_volume)

        btn_show_all = QPushButton("Show all profiles")
        btn_show_all.clicked.connect(self._on_show_all_profiles)
        left.addWidget(btn_show_all)

        btn_cam = QPushButton("Reset camera")
        btn_cam.clicked.connect(self._reset_camera)
        left.addWidget(btn_cam)

        if self.profiles:
            self.profile_list.setCurrentRow(0)

    def _on_wireframe_toggled(self, checked):
        """Slot called when Show wireframe checkbox changes."""
        self.show_wireframe = bool(checked)
        if hasattr(self, 'plotter'):
            if getattr(self, 'stl_file', None) is not None:
                self._show_stl()
            elif getattr(self, 'probes', None) is not None:
                self._show_probes()
            elif getattr(self, 'profiles', None) is not None:
                self._show_profiles()

    def _on_volume_toggled(self, checked):
        """Slot called when Translucent volume checkbox changes."""
        self.show_volume = bool(checked)
        if hasattr(self, 'plotter'):
            if getattr(self, 'stl_file', None) is not None:
                self._show_stl()
            elif getattr(self, 'probes', None) is not None:
                self._show_probes()
            elif getattr(self, 'profiles', None) is not None:
                self._show_profiles()

    def set_wireframe(self, enable):
        """Set wireframe display and re-render."""
        self.show_wireframe = bool(enable)
        if hasattr(self, 'check_wireframe'):
            self.check_wireframe.blockSignals(True)
            self.check_wireframe.setChecked(self.show_wireframe)
            self.check_wireframe.blockSignals(False)
        if hasattr(self, 'plotter'):
            if getattr(self, 'stl_file', None) is not None:
                self._show_stl()
            elif getattr(self, 'probes', None) is not None:
                self._show_probes()
            elif getattr(self, 'profiles', None) is not None:
                self._show_profiles()

    def set_volume(self, enable):
        """Set translucent volume display and re-render."""
        self.show_volume = bool(enable)
        if hasattr(self, 'check_volume'):
            self.check_volume.blockSignals(True)
            self.check_volume.setChecked(self.show_volume)
            self.check_volume.blockSignals(False)
        if hasattr(self, 'plotter'):
            if getattr(self, 'stl_file', None) is not None:
                self._show_stl()
            elif getattr(self, 'probes', None) is not None:
                self._show_probes()
            elif getattr(self, 'profiles', None) is not None:
                self._show_profiles()

    def _on_stl_toggled(self, checked):
        """Slot called when Show STL object checkbox changes."""
        self.show_stl = bool(checked)
        if hasattr(self, 'plotter'):
            self._show_stl()

    def _on_stl_rep_changed(self, index):
        """Slot called when STL representation combo changes."""
        if hasattr(self, 'combo_stl_rep'):
            self.stl_representation = self.combo_stl_rep.currentData() or "surface"
        if hasattr(self, 'plotter'):
            self._show_stl()

    def _on_seeds_toggled(self, checked):
        """Slot called when Show seed points checkbox changes."""
        self.show_seed_points = bool(checked)
        if hasattr(self, 'plotter'):
            self._show_stl()

    def _on_profile_display_mode_changed(self, index):
        """Slot called when profile display mode combo box changes."""
        mode = self.combo_profile_display.currentData()
        self.profile_display_mode = mode
        if hasattr(self, 'plotter'):
            self._show_profiles()

    def set_profile_display_mode(self, mode):
        """Set profile display mode ('line' or 'points') and re-render."""
        if mode not in ("line", "points", "spheres"):
            raise ValueError("Invalid display mode '%s'. Choose 'line' or 'points'." % mode)
        self.profile_display_mode = mode
        if hasattr(self, 'combo_profile_display'):
            idx = 1 if mode in ("points", "spheres") else 0
            self.combo_profile_display.blockSignals(True)
            self.combo_profile_display.setCurrentIndex(idx)
            self.combo_profile_display.blockSignals(False)
        if hasattr(self, 'plotter'):
            self._show_profiles()

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
            if hasattr(self, 'stl_polydata') and self.stl_polydata is not None:
                bounds = self.stl_polydata.bounds
            elif hasattr(self, 'grid') and self.grid is not None:
                bounds = self.grid.bounds
            else:
                return 0.01
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

    def _render_background_mesh(self):
        """
        Render boundary mesh background according to self.show_wireframe
        and self.show_volume settings.
        """
        if not (self.show_wireframe or self.show_volume):
            return
        if not hasattr(self, 'subsets') or not self.subsets:
            return

        for name, subset in self.subsets.items():
            if self.show_volume:
                self.plotter.add_mesh(subset,
                                       color="lightgray",
                                       opacity=0.15,
                                       show_edges=False,
                                       reset_camera=False)
            if self.show_wireframe:
                self.plotter.add_mesh(subset,
                                       style="wireframe",
                                       color="gray",
                                       reset_camera=False)

    def _refresh_probe_list(self):
        """Update QListWidget from self.probes (local list)."""
        self.probe_list.clear()
        for name, x, y, z in self.probes:
            self.probe_list.addItem(
                "%s  (%.3f, %.3f, %.3f)" % (name, x, y, z))

    def _show_probes(self):
        """
        Display mesh (wireframe and/or translucent volume) + all probes as spheres.

        Uses spheres precomputed in thread (_probe_spheres).
        Creates an individual actor per probe to support add/delete/highlight.
        """
        if not hasattr(self, 'plotter'):
            return

        cam = self._get_camera_state()
        self.plotter.clear()
        self._probe_actors       = {}
        self._selected_probe_idx = None

        # Boundary mesh background (wireframe and/or translucent volume)
        self._render_background_mesh()

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
    # Profiles mode
    # =========================================================================

    def _show_profiles(self, selected_idx=-1):
        """
        Display translucent mesh + all profiles as 3D curves/tubes inside volume.
        """
        if not hasattr(self, 'plotter'):
            return

        if selected_idx == -1:
            selected_idx = getattr(self, '_selected_profile_idx', None)

        cam = self._get_camera_state()
        self.plotter.clear()
        self._profile_actors = {}

        # Boundary mesh background (wireframe and/or translucent volume)
        self._render_background_mesh()

        r = self._probe_radius()
        tube_r = max(r * 0.35, 1e-6)
        sph_r  = max(r * 0.7, 1e-6)

        palette = [
            "#e74c3c",  # Red
            "#3498db",  # Blue
            "#2ecc71",  # Green
            "#f39c12",  # Orange
            "#9b59b6",  # Purple
            "#1abc9c",  # Turquoise
            "#e67e22",  # Carrot
            "#34495e",  # Slate
        ]

        mode = getattr(self, 'profile_display_mode', 'line')

        for idx, prof in enumerate(self.profiles):
            pts = np.asarray(prof.get('points', []), dtype=float)
            if len(pts) < 1:
                continue

            color = palette[idx % len(palette)]
            is_sel = (selected_idx is None or idx == selected_idx)
            if is_sel:
                tube_color  = color if selected_idx is None else "#e74c3c"
                start_color = "#2ecc71"  # Green for s=0
                end_color   = "#e74c3c"  # Red for s=1
                cur_tube_r  = tube_r * 1.4 if selected_idx is not None else tube_r
                cur_sph_r   = sph_r * 1.4 if selected_idx is not None else sph_r
            else:
                tube_color  = "#95a5a6"  # Subtle gray
                start_color = "#bdc3c7"
                end_color   = "#bdc3c7"
                cur_tube_r  = tube_r * 0.6
                cur_sph_r   = sph_r * 0.6

            try:
                if mode in ("points", "spheres"):
                    # Points mode: display one sphere per probe point defined in profile
                    acts = []
                    cur_point_r = cur_sph_r
                    if len(pts) > 1:
                        dists = np.linalg.norm(np.diff(pts, axis=0), axis=1)
                        min_dist = float(np.min(dists)) if len(dists) > 0 else 0.0
                        if min_dist > 1e-8:
                            cur_point_r = max(min(cur_sph_r, min_dist * 0.45), cur_sph_r * 0.25, 1e-6)

                    # Intermediate probe points
                    if len(pts) > 2:
                        mid_poly = pv.PolyData(pts[1:-1])
                        mid_glyph = mid_poly.glyph(geom=pv.Sphere(radius=cur_point_r), scale=False, orient=False)
                        act_mid = self.plotter.add_mesh(mid_glyph, color=tube_color, reset_camera=False)
                        acts.append(act_mid)
                    elif len(pts) == 1:
                        act_s0 = self.plotter.add_mesh(pv.Sphere(radius=cur_point_r, center=pts[0]),
                                                       color=start_color, reset_camera=False)
                        acts.append(act_s0)

                    if len(pts) >= 2:
                        act_s0 = self.plotter.add_mesh(pv.Sphere(radius=cur_point_r * 1.15, center=pts[0]),
                                                       color=start_color, reset_camera=False)
                        act_s1 = self.plotter.add_mesh(pv.Sphere(radius=cur_point_r * 1.15, center=pts[-1]),
                                                       color=end_color, reset_camera=False)
                        acts.extend([act_s0, act_s1])

                    self._profile_actors[idx] = tuple(acts)
                else:
                    if len(pts) < 2:
                        continue
                    line = pv.lines_from_points(pts)
                    tube = line.tube(radius=cur_tube_r)
                    act_tube = self.plotter.add_mesh(tube, color=tube_color, reset_camera=False)
                    act_s0   = self.plotter.add_mesh(pv.Sphere(radius=cur_sph_r, center=pts[0]),
                                                     color=start_color, reset_camera=False)
                    act_s1   = self.plotter.add_mesh(pv.Sphere(radius=cur_sph_r, center=pts[-1]),
                                                     color=end_color, reset_camera=False)
                    self._profile_actors[idx] = (act_tube, act_s0, act_s1)
            except Exception as e:
                log.warning("Could not render profile %s: %s", prof.get('label'), e)

        self._finalize_render(cam)

    def _on_profile_selected(self, row):
        """Highlight selected profile and update details panel."""
        if row < 0 or row >= len(self.profiles):
            return
        self._selected_profile_idx = row
        prof = self.profiles[row]
        pts = np.asarray(prof.get('points', []), dtype=float)

        info = "<b>%s</b><br>" % prof.get('label', 'Profile')
        if len(pts) >= 2:
            info += "<b>Start (s=0):</b> (%.4g, %.4g, %.4g)<br>" % (pts[0][0], pts[0][1], pts[0][2])
            info += "<b>End (s=1):</b> (%.4g, %.4g, %.4g)<br>" % (pts[-1][0], pts[-1][1], pts[-1][2])
            info += "<b>Points:</b> %d<br>" % len(pts)
        mode_str = "Points" if getattr(self, 'profile_display_mode', 'line') in ("points", "spheres") else "Straight line"
        info += "<b>Display:</b> %s<br>" % mode_str
        formula = prof.get('formula')
        if formula:
            clean_form = formula.strip().replace('\n', '; ')
            if len(clean_form) > 60:
                clean_form = clean_form[:57] + "..."
            info += "<b>Formula:</b> %s" % clean_form
        if hasattr(self, 'profile_details'):
            self.profile_details.setHtml(info)

        if hasattr(self, 'plotter'):
            self._show_profiles(selected_idx=row)

    def _on_show_all_profiles(self):
        """Reset profile selection to show all profiles uniformly."""
        self._selected_profile_idx = None
        if hasattr(self, 'profile_list'):
            self.profile_list.setCurrentRow(-1)
            self.profile_list.clearSelection()
        if hasattr(self, 'profile_details'):
            self.profile_details.setHtml("<i>All profiles displayed</i>")
        self._show_profiles(selected_idx=None)

    # =========================================================================
    # STL mode
    # =========================================================================

    def _show_stl(self):
        """
        Display background CFD mesh (wireframe and/or translucent volume) + STL object.
        """
        if not hasattr(self, 'plotter'):
            return

        cam = self._get_camera_state()
        self.plotter.clear()

        # Render background CFD mesh in wireframe and/or translucent volume
        self._render_background_mesh()

        # Load STL geometry if needed
        if self.stl_polydata is None and self.stl_file and os.path.isfile(self.stl_file):
            try:
                self.stl_polydata = pv.read(self.stl_file)
            except Exception as e:
                log.warning("Could not read STL file %s: %s", self.stl_file, e)

        # Render STL geometry
        if getattr(self, 'show_stl', True) and self.stl_polydata is not None:
            rep = getattr(self, 'stl_representation', 'surface')
            show_edges = False
            style = 'surface'
            if rep == 'surface_edges':
                show_edges = True
                style = 'surface'
            elif rep == 'wireframe':
                style = 'wireframe'
            elif rep == 'points':
                style = 'points'

            try:
                self.plotter.add_mesh(
                    self.stl_polydata,
                    color="#3498db",
                    style=style,
                    opacity=0.9,
                    show_edges=show_edges,
                    edge_color="black",
                    reset_camera=False
                )
            except Exception as e:
                log.warning("Error rendering STL mesh: %s", e)

        # Render seed points if any
        if getattr(self, 'show_seed_points', True) and self.seed_points:
            r = self._probe_radius()
            for item in self.seed_points:
                if len(item) == 4:
                    lbl, x, y, z = item
                elif len(item) == 3:
                    x, y, z = item
                    lbl = ""
                else:
                    continue
                try:
                    sph = pv.Sphere(radius=max(r, 1e-4), center=(float(x), float(y), float(z)))
                    self.plotter.add_mesh(sph, color="#e74c3c", reset_camera=False)
                except Exception as e:
                    log.warning("Error rendering seed point: %s", e)

        self._finalize_render(cam)

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
            if hasattr(self, 'profile_list'):
                self.profile_list.currentRowChanged.disconnect()
            if hasattr(self, 'combo_profile_display'):
                self.combo_profile_display.currentIndexChanged.disconnect()
            if hasattr(self, 'check_wireframe'):
                self.check_wireframe.toggled.disconnect()
            if hasattr(self, 'check_volume'):
                self.check_volume.toggled.disconnect()
        except Exception:
            pass  # disconnect might fail if signal was not connected

        # Wait for loading thread if still active
        if hasattr(self, '_loader') and self._loader.isRunning():
            self._loader.wait()

        # Cleanly close VTK plotter
        if hasattr(self, 'plotter'):
            self.plotter.close()

        event.accept()


#-------------------------------------------------------------------------------
# Unit tests
#-------------------------------------------------------------------------------

class MeshViewerDialogTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
        cls.app = QApplication.instance() or QApplication([])

    def testAxesWidgetConfiguration(self):
        """
        Verify that QtInteractor in MeshViewerDialog configures orientation axes
        (X, Y, Z) in the lower left corner.
        """
        if not HAS_VISU:
            self.skipTest("PyVista/PyVistaQt not available")

        dlg = MeshViewerDialog.__new__(MeshViewerDialog)
        super(MeshViewerDialog, dlg).__init__(None)
        dlg.probes = None
        dlg.profiles = None
        dlg.name_to_tags = {'Wall': [1]}
        dlg.zone_expr = None
        dlg.zone_label = None
        dlg.case = None
        dlg._probes_modified = False
        dlg._unique_names = ['Wall']

        dlg._build_ui()

        # Check that axes widget is present and enabled
        widget = getattr(dlg.plotter.renderer, 'axes_widget', None)
        self.assertIsNotNone(widget)
        self.assertEqual(widget.GetEnabled(), 1)

        # Check lower-left viewport: (0, 0, xmax, ymax) with positive dimensions
        vp = widget.GetViewport()
        self.assertAlmostEqual(vp[0], 0.0)
        self.assertAlmostEqual(vp[1], 0.0)
        self.assertGreater(vp[2], 0.0)
        self.assertGreater(vp[3], 0.0)

        # Check axes labels: X, Y, Z
        axes_actor = widget.GetOrientationMarker()
        self.assertIsNotNone(axes_actor)
        self.assertEqual(axes_actor.GetXAxisLabelText(), 'X')
        self.assertEqual(axes_actor.GetYAxisLabelText(), 'Y')
        self.assertEqual(axes_actor.GetZAxisLabelText(), 'Z')

        dlg.close()

    def testEvalProfileFormula(self):
        """
        Verify evaluation of profile mathematical formulas with parameter s in [0, 1].
        """
        formula = """
        x0 = 1.0; x1 = 4.0;
        y0 = 2.0; y1 = 6.0;
        z0 = 3.0; z1 = 8.0;
        x = x0 + (x1 - x0)*s;
        y = y0 + (y1 - y0)*s;
        z = z0 + (z1 - z0)*s;
        """
        pts = eval_profile_formula(formula, nb_points=10)
        self.assertEqual(len(pts), 10)
        np.testing.assert_allclose(pts[0], [1.0, 2.0, 3.0])
        np.testing.assert_allclose(pts[-1], [4.0, 6.0, 8.0])

        # Test formula with notebook variables
        f_nb = "x = L * s; y = 0.0; z = 0.0;"
        pts_nb = eval_profile_formula(f_nb, nb_points=5, notebook_vars={'L': 12.5})
        self.assertEqual(len(pts_nb), 5)
        np.testing.assert_allclose(pts_nb[0], [0.0, 0.0, 0.0])
        np.testing.assert_allclose(pts_nb[-1], [12.5, 0.0, 0.0])

        # Test invalid formula raises ValueError
        with self.assertRaises(ValueError):
            eval_profile_formula("x = s; y = s;")  # missing z

    def testProfileModeUI(self):
        """
        Verify MeshViewerDialog UI initialization in profile visualization mode.
        """
        if not HAS_VISU:
            self.skipTest("PyVista/PyVistaQt not available")

        pts = np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]], dtype=float)
        prof = {'label': 'Profile_1', 'points': pts, 'formula': 'x=s; y=s; z=s;', 'nb_points': 2}

        dlg = MeshViewerDialog.__new__(MeshViewerDialog)
        super(MeshViewerDialog, dlg).__init__(None)
        dlg.probes = None
        dlg.profiles = [prof]
        dlg.name_to_tags = {'Wall': [1]}
        dlg.zone_expr = None
        dlg.zone_label = None
        dlg.case = None
        dlg._probes_modified = False
        dlg._unique_names = ['Wall']

        dlg._build_ui()

        self.assertTrue(hasattr(dlg, 'profile_list'))
        self.assertEqual(dlg.profile_list.count(), 1)
        self.assertEqual(dlg.profile_list.item(0).text(), 'Profile_1')

        # Test rendering profiles
        dlg.subsets = {}
        dlg._show_profiles()
        self.assertIn(0, dlg._profile_actors)

        # Test profile details and selection
        self.assertIn('Profile_1', dlg.profile_details.toHtml())
        dlg._on_show_all_profiles()
        self.assertIsNone(dlg._selected_profile_idx)
        self.assertIn('All profiles displayed', dlg.profile_details.toHtml())
        dlg._on_profile_selected(0)
        self.assertEqual(dlg._selected_profile_idx, 0)
        self.assertIn('Profile_1', dlg.profile_details.toHtml())

        dlg.close()

    def testProfileModePointsUI(self):
        """
        Verify MeshViewerDialog points representation mode for profiles.
        """
        if not HAS_VISU:
            self.skipTest("PyVista/PyVistaQt not available")

        # 5 discrete probe points along a segment
        pts = np.array([
            [0.0, 0.0, 0.0],
            [0.25, 0.25, 0.25],
            [0.5, 0.5, 0.5],
            [0.75, 0.75, 0.75],
            [1.0, 1.0, 1.0]
        ], dtype=float)
        prof = {'label': 'Profile_Points', 'points': pts, 'formula': 'x=s; y=s; z=s;', 'nb_points': 5}

        dlg = MeshViewerDialog.__new__(MeshViewerDialog)
        super(MeshViewerDialog, dlg).__init__(None)
        dlg.probes = None
        dlg.profiles = [prof]
        dlg.profile_display_mode = "points"
        dlg.name_to_tags = {'Wall': [1]}
        dlg.zone_expr = None
        dlg.zone_label = None
        dlg.case = None
        dlg._probes_modified = False
        dlg._unique_names = ['Wall']

        dlg._build_ui()

        # Check UI combo box exists and points mode is selected
        self.assertTrue(hasattr(dlg, 'combo_profile_display'))
        self.assertEqual(dlg.combo_profile_display.currentData(), 'points')
        self.assertEqual(dlg.combo_profile_display.currentIndex(), 1)

        # Test rendering profiles in points mode
        dlg.subsets = {}
        dlg._show_profiles()
        self.assertIn(0, dlg._profile_actors)
        # Should have intermediate glyph actor + start + end sphere actors
        actors = dlg._profile_actors[0]
        self.assertGreaterEqual(len(actors), 2)

        # Check profile details contains Points mode indicator
        dlg._on_profile_selected(0)
        self.assertIn('Points', dlg.profile_details.toHtml())

        # Test switching to line mode
        dlg.set_profile_display_mode("line")
        self.assertEqual(dlg.profile_display_mode, "line")
        self.assertEqual(dlg.combo_profile_display.currentData(), "line")
        self.assertEqual(dlg.combo_profile_display.currentIndex(), 0)
        self.assertIn(0, dlg._profile_actors)
        line_actors = dlg._profile_actors[0]
        self.assertEqual(len(line_actors), 3)  # tube, s0, s1

        # Test switching back to points mode
        dlg.set_profile_display_mode("points")
        self.assertEqual(dlg.profile_display_mode, "points")
        self.assertEqual(dlg.combo_profile_display.currentData(), "points")

        # Test invalid mode raises ValueError
        with self.assertRaises(ValueError):
            dlg.set_profile_display_mode("invalid_mode")

        dlg.close()

    def testProbeModeWireframeUI(self):
        """
        Verify MeshViewerDialog UI and rendering controls in probe mode with wireframe.
        """
        if not HAS_VISU:
            self.skipTest("PyVista/PyVistaQt not available")

        dlg = MeshViewerDialog.__new__(MeshViewerDialog)
        super(MeshViewerDialog, dlg).__init__(None)
        dlg.probes = [('probe_1', 0.1, 0.2, 0.3)]
        dlg.profiles = None
        dlg.show_wireframe = True
        dlg.show_volume = True
        dlg.name_to_tags = {'Wall': [1]}
        dlg.zone_expr = None
        dlg.zone_label = None
        dlg.case = None
        dlg._probes_modified = False
        dlg._unique_names = ['Wall']

        dlg._build_ui()

        # Check UI display options checkboxes
        self.assertTrue(hasattr(dlg, 'check_wireframe'))
        self.assertTrue(hasattr(dlg, 'check_volume'))
        self.assertTrue(dlg.check_wireframe.isChecked())
        self.assertTrue(dlg.check_volume.isChecked())

        # Test rendering in probe mode
        dlg.subsets = {}
        dlg._probe_spheres = []
        dlg._show_probes()
        self.assertIn(0, dlg._probe_actors)

        # Toggle wireframe off
        dlg.check_wireframe.setChecked(False)
        self.assertFalse(dlg.show_wireframe)

        # Toggle volume off
        dlg.check_volume.setChecked(False)
        self.assertFalse(dlg.show_volume)

        # Programmatic setters
        dlg.set_wireframe(True)
        self.assertTrue(dlg.show_wireframe)
        self.assertTrue(dlg.check_wireframe.isChecked())

        dlg.set_volume(True)
        self.assertTrue(dlg.show_volume)
        self.assertTrue(dlg.check_volume.isChecked())

        dlg.close()

    def testProfileModeWireframeUI(self):
        """
        Verify MeshViewerDialog UI and rendering controls in profile mode with wireframe.
        """
        if not HAS_VISU:
            self.skipTest("PyVista/PyVistaQt not available")

        pts = np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]], dtype=float)
        prof = {'label': 'Profile_Wireframe', 'points': pts, 'formula': 'x=s; y=s; z=s;', 'nb_points': 2}

        dlg = MeshViewerDialog.__new__(MeshViewerDialog)
        super(MeshViewerDialog, dlg).__init__(None)
        dlg.probes = None
        dlg.profiles = [prof]
        dlg.show_wireframe = True
        dlg.show_volume = True
        dlg.name_to_tags = {'Wall': [1]}
        dlg.zone_expr = None
        dlg.zone_label = None
        dlg.case = None
        dlg._probes_modified = False
        dlg._unique_names = ['Wall']

        dlg._build_ui()

        # Check UI display options checkboxes
        self.assertTrue(hasattr(dlg, 'check_wireframe'))
        self.assertTrue(hasattr(dlg, 'check_volume'))
        self.assertTrue(dlg.check_wireframe.isChecked())
        self.assertTrue(dlg.check_volume.isChecked())

        # Test rendering in profile mode
        dlg.subsets = {}
        dlg._show_profiles()
        self.assertIn(0, dlg._profile_actors)

        # Toggle wireframe off
        dlg.check_wireframe.setChecked(False)
        self.assertFalse(dlg.show_wireframe)

        # Toggle volume off
        dlg.check_volume.setChecked(False)
        self.assertFalse(dlg.show_volume)

        # Programmatic setters
        dlg.set_wireframe(True)
        self.assertTrue(dlg.show_wireframe)
        self.assertTrue(dlg.check_wireframe.isChecked())

        dlg.set_volume(True)
        self.assertTrue(dlg.show_volume)
        self.assertTrue(dlg.check_volume.isChecked())

        dlg.close()

    def testSTLModeUI(self):
        """
        Verify MeshViewerDialog UI controls in STL mode.
        """
        if not HAS_VISU:
            self.skipTest("PyVista/PyVistaQt not available")

        import tempfile
        cube = pv.Cube()
        with tempfile.NamedTemporaryFile(suffix='.stl', delete=False) as f:
            stl_path = f.name
        try:
            cube.save(stl_path)

            dlg = MeshViewerDialog.__new__(MeshViewerDialog)
            super(MeshViewerDialog, dlg).__init__(None)
            dlg.csm_file = None
            dlg.probes = None
            dlg.profiles = None
            dlg.stl_file = stl_path
            dlg.stl_label = "TestCube"
            dlg.seed_points = [("Seed_1", 0.0, 0.0, 1.0)]
            dlg.show_wireframe = True
            dlg.show_volume = True
            dlg.show_stl = True
            dlg.show_stl_edges = False
            dlg.show_seed_points = True
            dlg.stl_representation = "surface"
            dlg.name_to_tags = {}
            dlg.zone_expr = None
            dlg.zone_label = None
            dlg.case = None
            dlg._probes_modified = False
            dlg._unique_names = []

            dlg._build_ui()

            self.assertTrue(hasattr(dlg, 'check_stl'))
            self.assertTrue(hasattr(dlg, 'combo_stl_rep'))
            self.assertTrue(hasattr(dlg, 'check_wireframe'))
            self.assertTrue(hasattr(dlg, 'check_volume'))
            self.assertTrue(hasattr(dlg, 'check_seeds'))
            self.assertTrue(hasattr(dlg, 'stl_details'))

            self.assertTrue(dlg.check_stl.isChecked())
            self.assertTrue(dlg.check_wireframe.isChecked())
            self.assertTrue(dlg.check_volume.isChecked())
            self.assertTrue(dlg.check_seeds.isChecked())

            # Test toggling STL visibility
            dlg.check_stl.setChecked(False)
            self.assertFalse(dlg.show_stl)
            dlg.check_stl.setChecked(True)
            self.assertTrue(dlg.show_stl)

            # Test toggling representation
            dlg.combo_stl_rep.setCurrentIndex(1)
            self.assertEqual(dlg.stl_representation, "surface_edges")

            dlg.close()
        finally:
            if os.path.exists(stl_path):
                os.remove(stl_path)

    def testSTLPositionCheck(self):
        """
        Verify automated bounds check: inside vs outside mesh bounds.
        """
        if not HAS_VISU:
            self.skipTest("PyVista/PyVistaQt not available")

        import tempfile
        # Cube of size 1 centered at origin: bounds [-0.5, 0.5, -0.5, 0.5, -0.5, 0.5]
        cube = pv.Cube()
        with tempfile.NamedTemporaryFile(suffix='.stl', delete=False) as f:
            stl_path = f.name
        try:
            cube.save(stl_path)

            dlg = MeshViewerDialog.__new__(MeshViewerDialog)
            super(MeshViewerDialog, dlg).__init__(None)
            dlg.csm_file = "mesh.csm"
            dlg.probes = None
            dlg.profiles = None
            dlg.stl_file = stl_path
            dlg.stl_label = "InsideCube"
            dlg.seed_points = []
            dlg.show_wireframe = True
            dlg.show_volume = True
            dlg.show_stl = True
            dlg.show_stl_edges = False
            dlg.show_seed_points = False
            dlg.stl_representation = "surface"
            dlg.name_to_tags = {'Domain': [1]}
            dlg.zone_expr = None
            dlg.zone_label = None
            dlg.case = None
            dlg._probes_modified = False
            dlg._unique_names = ['Domain']

            # Background domain mesh encompassing the cube: [-2.0, 2.0, -2.0, 2.0, -2.0, 2.0]
            domain_cube = pv.Cube(bounds=[-2.0, 2.0, -2.0, 2.0, -2.0, 2.0])
            dlg.grid = domain_cube
            dlg.subsets = {'Domain': domain_cube}

            dlg._build_ui()

            self.assertTrue(hasattr(dlg, 'label_position_check'))
            self.assertIn("INSIDE", dlg.label_position_check.text())
            dlg.close()

            # Now test an STL extending outside the domain: domain is [-0.1, 0.1, ...]
            dlg2 = MeshViewerDialog.__new__(MeshViewerDialog)
            super(MeshViewerDialog, dlg2).__init__(None)
            dlg2.csm_file = "mesh.csm"
            dlg2.probes = None
            dlg2.profiles = None
            dlg2.stl_file = stl_path
            dlg2.stl_label = "OutsideCube"
            dlg2.seed_points = []
            dlg2.show_wireframe = True
            dlg2.show_volume = True
            dlg2.show_stl = True
            dlg2.show_stl_edges = False
            dlg2.show_seed_points = False
            dlg2.stl_representation = "surface"
            dlg2.name_to_tags = {'Domain': [1]}
            dlg2.zone_expr = None
            dlg2.zone_label = None
            dlg2.case = None
            dlg2._probes_modified = False
            dlg2._unique_names = ['Domain']

            small_domain = pv.Cube(bounds=[-0.1, 0.1, -0.1, 0.1, -0.1, 0.1])
            dlg2.grid = small_domain
            dlg2.subsets = {'Domain': small_domain}

            dlg2._build_ui()
            self.assertTrue(hasattr(dlg2, 'label_position_check'))
            self.assertIn("OUTSIDE", dlg2.label_position_check.text())
            dlg2.close()
        finally:
            if os.path.exists(stl_path):
                os.remove(stl_path)


def suite():
    testSuite = unittest.TestSuite()
    testSuite.addTest(unittest.TestLoader().loadTestsFromTestCase(MeshViewerDialogTestCase))
    return testSuite


def runTest():
    runner = unittest.TextTestRunner()
    runner.run(suite())


if __name__ == '__main__':
    runTest()
