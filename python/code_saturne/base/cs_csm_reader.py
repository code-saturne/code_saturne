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
This module provides functions and classes to read Code_Saturne native
mesh files (.csm format, Big Endian binary).

.csm file structure (Big Endian):

  Global header (152 bytes):
    [0-63]   : magic base  "Code_Saturne I/O, BE, R0" (padded to 64)
    [64-127] : magic file  "Face-based mesh definition, R0" (padded to 64)
    [128-135]: header_size  uint64 BE — fixed size of each section header
                            (96 or 128 depending on Code_Saturne version)
    [136-143]: header_align uint64 BE — header alignment (typically 64)
    [144-151]: body_align   uint64 BE — data alignment (typically 64)

  Section header (header_size bytes):
    [0-7]  : total_size — actual header size (may exceed header_size
                          if section name is long)
    [8-15] : n_vals     — number of values in this section
    [16-23]: location_id (unused here)
    [24-31]: index_id    (unused here)
    [32-39]: n_loc_vals  (unused here)
    [48-49]: type string — 'r4','r8','i4','i8','u4','u8','c ' or '  '
    [55]   : flag 'e' if data is embedded in the header
    [56+]  : section name (null-terminated, aligned to 8 bytes)

  Key sections for visualization:
    face_cells          : (n_faces, 2) — cell[1]==0 → boundary face
    face_group_class_id : (n_faces,)   — group class of each face
    face_vertices_index : (n_faces+1,) — CSR index of vertices
    face_vertices       : (n_verts,)   — face→vertices connectivity (1-based)
    vertex_coords       : (n_pts×3,)   — XYZ coordinates
    group_name          : string       — group names separated by \x00
    group_class_properties : (n_props×n_gc,) — gc_id → name mapping
"""

import os
import re
import struct
import unittest
import numpy as np

# Byte size of each data type
TYPE_SIZES = {
    'r4': 4, 'r8': 8,
    'i4': 4, 'i8': 8,
    'u4': 4, 'u8': 8,
    'c ': 1,           # string
    '  ': 0,           # control section (no data)
}

# CSM type to numpy dtype mapping (Big Endian)
DTYPE_MAP = {
    'r4': '>f4', 'r8': '>f8',
    'i4': '>i4', 'i8': '>i8',
    'u4': '>u4', 'u8': '>u8',
}


def _align(pos, alignment):
    """Return the next position aligned to `alignment` bytes."""
    if alignment == 0:
        return pos
    return pos + (alignment - pos % alignment) % alignment


def read_csm(csm_file):
    """
    Read a Code_Saturne .csm file (Big Endian binary format).

    Returns a dictionary {section_name: data} where data is either
    a np.ndarray (for numeric types) or a str (for 'c ').

    Handles both known variants of the format:
    - header_size = 128 (recent Code_Saturne versions)
    - header_size = 96  (older versions)

    Special cases handled:
    - total_size > header_size: long section name → reads surplus
    - total_size < header_size: control section → no negative seek
    """
    result = {}

    with open(csm_file, 'rb') as f:

        #  Global header reading
        f.read(64)  # magic base (ignored, file extension trusted)
        f.read(64)  # magic file (ignored)
        header_size  = struct.unpack('>Q', f.read(8))[0]
        header_align = struct.unpack('>Q', f.read(8))[0]
        body_align   = struct.unpack('>Q', f.read(8))[0]

        #  Loop over sections
        while True:

            # Alignment before each section header
            pos = f.tell()
            new_pos = _align(pos, header_align)
            if new_pos > pos:
                f.seek(new_pos)

            # Read section header
            header_raw = f.read(header_size)
            if len(header_raw) < header_size:
                break  # EOF

            # Actual header size (can exceed header_size if name is long)
            total_size = struct.unpack_from('>Q', header_raw, 0)[0]

            # If actual header is larger, read the surplus
            # (e.g. long section name like "group_name" with 132 chars)
            if 0 < total_size <= 65536 and total_size > header_size:
                extra_raw  = f.read(total_size - header_size)
                header_raw = header_raw + extra_raw

            # Decode header
            n_vals   = struct.unpack_from('>Q', header_raw,  8)[0]
            type_str = header_raw[48:50].decode('ascii', errors='replace')
            # Flag 'e' at byte 55: data is embedded in header
            embedded = (header_raw[55] == ord('e'))

            # Section name (null-terminated starting at byte 56)
            name_raw = header_raw[56:]
            null_pos = name_raw.find(b'\x00')
            sec_name = (name_raw[:null_pos].decode('ascii', errors='replace')
                        if null_pos >= 0
                        else name_raw.decode('ascii', errors='replace'))

            if sec_name == 'EOF':
                break

            type_size = TYPE_SIZES.get(type_str, 0)

            # Control section (start_block, end_block, etc.) → skip
            # Important: do not seek, header_size bytes have already been read
            if n_vals == 0 or type_size == 0:
                continue

            #  Data reading
            if embedded:
                # Data in the header itself, after the null-terminated name
                # aligned to 8 bytes
                name_end = 56
                while name_end < len(header_raw) and header_raw[name_end] != 0:
                    name_end += 1
                name_end += 1                    # skip \x00
                name_end = (name_end + 7) & ~7   # align to 8 bytes
                data_raw = header_raw[name_end:name_end + n_vals * type_size]
            else:
                # Data following header, aligned to body_align
                pos = f.tell()
                new_pos = _align(pos, body_align)
                if new_pos > pos:
                    f.seek(new_pos)
                data_raw = f.read(n_vals * type_size)

            #  Storage
            if type_str == 'c ':
                # ASCII string
                result[sec_name] = (data_raw
                                    .decode('ascii', errors='replace')
                                    .rstrip('\x00'))
            elif type_str in DTYPE_MAP:
                # Numpy array (copied to avoid buffer ownership issues)
                result[sec_name] = np.frombuffer(
                    data_raw, dtype=DTYPE_MAP[type_str]).copy()

    return result


def decode_csm_groups(data):
    """
    Decode group names and group class mappings from CSM file data.

    Returns
    -------
    gc_to_name : dict {int: str}
        Mapping gc_id (1-based) → group name
    name_to_tags : dict {str: list[int]}
        Mapping group name → list of gc_ids
    """
    group_name_raw = data.get('group_name', '')
    group_names = [g for g in group_name_raw.split('\x00') if g]

    n_gc = int(data['n_group_classes'][0]) if 'n_group_classes' in data else 0
    n_props_max = int(data['n_group_class_props_max'][0]) if 'n_group_class_props_max' in data else 0

    gc_props_flat = data.get('group_class_properties', np.array([], dtype=int))
    if n_props_max > 0 and len(gc_props_flat) % n_props_max == 0:
        n_gc_real = len(gc_props_flat) // n_props_max
    else:
        n_gc_real = n_gc

    if n_props_max > 0 and n_gc_real > 0:
        gc_props = gc_props_flat.reshape(n_props_max, n_gc_real)
    else:
        gc_props = np.zeros((0, 0), dtype=int)
    n_gc = n_gc_real

    gc_to_name = {}
    name_to_tags = {}

    for gc_idx in range(n_gc):
        gc_id = gc_idx + 1  # gc_id is 1-based
        names = []
        for prop_idx in range(n_props_max):
            val = gc_props[prop_idx, gc_idx]
            if val < 0:
                name_idx = (-val) - 1  # 1-based → 0-based
                if name_idx < len(group_names):
                    names.append(group_names[name_idx])
        # gc_id without name → default family (unassigned faces)
        name = " + ".join(names) if names else "Default (untagged)"
        gc_to_name[gc_id] = name
        if name not in name_to_tags:
            name_to_tags[name] = []
        name_to_tags[name].append(gc_id)

    return gc_to_name, name_to_tags


def csm_to_pyvista_boundary(csm_file):
    """
    Read a .csm file and build pyvista mesh of boundary faces.

    Boundary face identification:
        A face is a boundary face if face_cells[i][1] == 0
        (the face is adjacent to only one cell).

    Returns
    -------
    grid : pv.UnstructuredGrid
        Boundary face mesh with cell_data["tag"] = gc_id
    gc_to_name : dict {int: str}
        Mapping gc_id → group name
    name_to_tags : dict {str: list[int]}
        Mapping name → list of gc_ids (for composite groups)
    subsets : dict {str: pv.UnstructuredGrid}
        Pre-extracted submesh per group (used in _show_all)
    group_pairs : dict {str: (pv.UnstructuredGrid, pv.UnstructuredGrid)}
        Pair (group_faces, rest) pre-extracted per group
        (used in _on_group_selected → max 2 add_mesh calls on click)
    """
    try:
        import pyvista as pv
    except ImportError:
        raise ImportError(
            "The 'pyvista' library is required to generate the visualization grid."
        )

    data = read_csm(csm_file)

    #  Vertex coordinates
    points = data['vertex_coords'].astype(np.float64).reshape(-1, 3)

    #  Boundary face identification
    # face_cells is an array (n_faces, 2) where cell[1]==0 → boundary face
    face_cells  = data['face_cells'].reshape(-1, 2)
    face_gc_id  = data['face_group_class_id']

    # Face→vertices connectivity (1-based indices in CSM → 0-based)
    fv_index    = data['face_vertices_index'].astype(np.int64) - 1
    fv          = data['face_vertices'].astype(np.int64) - 1

    #  Group names decoding
    gc_to_name, name_to_tags = decode_csm_groups(data)

    #  grid construction
    boundary_mask     = face_cells[:, 1] == 0
    boundary_face_ids = np.where(boundary_mask)[0]
    boundary_gc       = face_gc_id[boundary_face_ids]

    cells_list     = []
    celltypes_list = []
    tags_list      = []

    for i, face_id in enumerate(boundary_face_ids):
        start = fv_index[face_id]
        end   = fv_index[face_id + 1]
        verts = fv[start:end]
        n     = len(verts)

        # Determine pyvista cell type based on vertex count
        if n == 3:
            ct = pv.CellType.TRIANGLE
        elif n == 4:
            ct = pv.CellType.QUAD
        else:
            ct = pv.CellType.POLYGON  # general polygon (polyhedron faces)

        # format: [n_pts, p0, p1, ..., pn]
        cells_list.append(np.concatenate([[n], verts]))
        celltypes_list.append(ct)
        tags_list.append(boundary_gc[i])

    ct_np   = np.array(celltypes_list)
    tags_np = np.array(tags_list)

    # Remap points: keep only points referenced by faces
    # (.csm contains all volume points, ~10x more than needed)
    all_used        = np.unique(np.concatenate([c[1:] for c in cells_list]))
    remap           = np.zeros(len(points), dtype=np.int64)
    remap[all_used] = np.arange(len(all_used))
    points_small    = points[all_used]

    new_cells = [np.concatenate([[c[0]], remap[c[1:]]]) for c in cells_list]
    cells_np  = np.concatenate(new_cells)

    grid = pv.UnstructuredGrid(cells_np, ct_np, points_small)
    grid.cell_data["tag"] = tags_np

    # Subset pre-extraction (rendering optimization)
    # These extractions are performed in the loading thread to avoid
    # any computation in the main Qt thread during user interactions.

    # subsets: one grid per group → used in _show_all (N add_mesh calls)
    subsets = {}
    for name, tag_list in name_to_tags.items():
        mask = np.isin(grid.cell_data["tag"], tag_list)
        subset = grid.extract_cells(mask)
        if subset.n_cells > 0:
            subsets[name] = subset

    # group_pairs: (selected, rest) per group → used in _on_group_selected
    # Allows calling add_mesh only twice instead of N times on highlight
    group_pairs = {}
    for name, tag_list in name_to_tags.items():
        mask_sel  = np.isin(grid.cell_data["tag"], tag_list)
        sel_grid  = grid.extract_cells(mask_sel)
        rest_grid = grid.extract_cells(~mask_sel)
        if sel_grid.n_cells > 0:
            group_pairs[name] = (sel_grid, rest_grid)

    return grid, gc_to_name, name_to_tags, subsets, group_pairs


def parse_localization(expr, grid, name_to_tags):
    """
    Parse a boundary zone expression and return a boolean numpy mask.

    Supported syntax:
        - Group name: "WALL1"
        - OR operator:  "WALL1 or WALL2"
        - AND operator: "WALL1 and WALL2"
        - NOT operator: "not WALL1" or "not(WALL1 or WALL2)"
        - Parentheses:  "(WALL1 or WALL2) and not(WALL3)"
        - Select all:   "all[]"
    """
    tags    = grid.cell_data["tag"]
    n_cells = len(tags)

    if expr.strip().lower() == "all[]":
        return np.ones(n_cells, dtype=bool)

    def get_mask_for_name(name):
        gc_ids = name_to_tags.get(name)
        if gc_ids is not None:
            return np.isin(tags, gc_ids)
        return np.zeros(n_cells, dtype=bool)

    def eval_expr(tokens, pos):
        """
        Recursively evaluate an expression starting from pos.
        Returns (mask, new_position).
        """
        mask = None
        op   = "or"

        while pos < len(tokens):
            tok = tokens[pos]

            if tok == ')':
                # End of parenthesized group
                return mask, pos + 1

            elif tok.lower() == 'or':
                op = "or"
                pos += 1

            elif tok.lower() == 'and':
                op = "and"
                pos += 1

            elif tok.lower() == 'not':
                pos += 1
                if pos < len(tokens) and tokens[pos] == '(':
                    # not(...): evaluate subgroup and negate result
                    pos += 1  # skip '('
                    sub_mask, pos = eval_expr(tokens, pos)
                    current = ~sub_mask if sub_mask is not None \
                              else np.ones(n_cells, dtype=bool)
                elif pos < len(tokens):
                    # not TOKEN: negate following token
                    current = ~get_mask_for_name(tokens[pos])
                    pos += 1
                else:
                    continue

                mask = current if mask is None else (
                    mask | current if op == "or" else mask & current)

            elif tok == '(':
                # Parenthesized group
                pos += 1
                current, pos = eval_expr(tokens, pos)
                if current is not None:
                    mask = current if mask is None else (
                        mask | current if op == "or" else mask & current)

            else:
                # Group name
                current = get_mask_for_name(tok)
                mask = current if mask is None else (
                    mask | current if op == "or" else mask & current)
                pos += 1

        return mask, pos

    # Tokenization: also split on parentheses
    raw_tokens = re.split(r'\s+', expr.strip())
    tokens = []
    for t in raw_tokens:
        # Handle cases like "not(WALL1" or "(WALL1" without space
        t = t.replace('(', ' ( ').replace(')', ' ) ')
        tokens.extend(t.split())
    tokens = [t for t in tokens if t]  # remove empty tokens

    mask, _ = eval_expr(tokens, 0)
    return mask if mask is not None else np.zeros(n_cells, dtype=bool)


class CsmReader(object):
    """
    Helper class to read and inspect a Code_Saturne .csm mesh file.
    """

    def __init__(self, csm_file):
        """
        Constructor.
        """
        self.csm_file = csm_file
        self._data = None
        self._groups = None

    def read(self):
        """
        Read raw CSM data dictionary.
        """
        if self._data is None:
            self._data = read_csm(self.csm_file)
        return self._data

    def get_groups(self):
        """
        Get (gc_to_name, name_to_tags) group dictionaries.
        """
        if self._groups is None:
            self._groups = decode_csm_groups(self.read())
        return self._groups

    def to_pyvista_boundary(self):
        """
        Build pyvista boundary grid and pre-extracted subsets.
        """
        return csm_to_pyvista_boundary(self.csm_file)


#-------------------------------------------------------------------------------
# Unit tests
#-------------------------------------------------------------------------------

class CsmReaderTestCase(unittest.TestCase):
    def testAlign(self):
        self.assertEqual(_align(0, 64), 0)
        self.assertEqual(_align(1, 64), 64)
        self.assertEqual(_align(64, 64), 64)
        self.assertEqual(_align(65, 64), 128)
        self.assertEqual(_align(10, 0), 10)

    def testDecodeGroups(self):
        data = {
            'group_name': 'Inlet\x00Outlet\x00Wall\x00',
            'n_group_classes': np.array([3], dtype=np.int32),
            'n_group_class_props_max': np.array([1], dtype=np.int32),
            'group_class_properties': np.array([-1, -2, -3], dtype=np.int32),
        }
        gc_to_name, name_to_tags = decode_csm_groups(data)
        self.assertEqual(gc_to_name[1], 'Inlet')
        self.assertEqual(gc_to_name[2], 'Outlet')
        self.assertEqual(gc_to_name[3], 'Wall')
        self.assertEqual(name_to_tags['Inlet'], [1])
        self.assertEqual(name_to_tags['Outlet'], [2])
        self.assertEqual(name_to_tags['Wall'], [3])

    def testParseLocalization(self):
        class MockGrid:
            def __init__(self):
                self.cell_data = {'tag': np.array([1, 2, 3, 1, 2])}
        grid = MockGrid()
        name_to_tags = {'Inlet': [1], 'Outlet': [2], 'Wall': [3]}

        m_all = parse_localization('all[]', grid, name_to_tags)
        self.assertTrue(m_all.all())

        m_or = parse_localization('Inlet or Outlet', grid, name_to_tags)
        np.testing.assert_array_equal(m_or, np.array([True, True, False, True, True]))

        m_not = parse_localization('not Wall', grid, name_to_tags)
        np.testing.assert_array_equal(m_not, np.array([True, True, False, True, True]))

        m_and = parse_localization('Inlet and Wall', grid, name_to_tags)
        self.assertFalse(m_and.any())


def suite():
    testSuite = unittest.TestSuite()
    testSuite.addTest(unittest.TestLoader().loadTestsFromTestCase(CsmReaderTestCase))
    return testSuite


def runTest():
    runner = unittest.TextTestRunner()
    runner.run(suite())


if __name__ == '__main__':
    runTest()
