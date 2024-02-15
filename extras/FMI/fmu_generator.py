#!/usr/bin/python3
# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2024 EDF S.A.
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

import os
import sys

real_vars = []
integer_vars = []
bool_vars = []
string_vars = []

def read_from_xml(path):

    import xml.etree.ElementTree as ET
    tree = ET.parse(path)
    root = tree.getroot()

    idx = 0
    for var in root.iter('ScalarVariable'):
        # print(var.attrib)
        d = {}
        d['index'] = idx
        for k in var.attrib.keys():
            d[k] = var.attrib[k]
        for child in var:
            for k in child.attrib.keys():
                d[k] = child.attrib[k]
            if child.tag == "Real":
                real_vars.append(d)
            elif child.tag == "Int":
                integer_vars.append(d)
            elif child.tag == "Int":
                bool_vars.append(d)
            elif child.tag == "String":
                string_vars.append(d)
        idx += 1

#-------------------------------------------------------------------------------

def generate_include_file(path):
    """
    Generate C/C++ file to include in FMU
    """

    lines = []

    lines.append("#pragma once")
    lines.append('')

    n_real = len(real_vars)
    n_int = len(integer_vars)
    n_bool = len(bool_vars)
    n_string = len(string_vars)

    lines.append("static const int _n_real_vars = {};".format(n_real))
    lines.append("static const int _n_integer_vars = {};".format(n_int))
    lines.append("static const int _n_bool_vars = {};".format(n_bool))
    lines.append("static const int _n_string_vars = {};".format(n_string))
    lines.append('')

    # Not all names currently needed
    # for var_type in ("real", "integer", "bool", "string"):
    for var_type in ("real", "string"):
        lines.append('static const char * _' + var_type + '_names[] = {')
        for var in eval(var_type + '_vars'):
            lines.append('  "' + var['name'] + '",')
        lines.append('};')
        lines.append('')

    for var_type in ("real", "integer", "bool"):
        if var_type == 'real':
            n_vars = n_real
            start_undef = '-HUGE_VAL_F'
        elif var_type == 'integer':
            n_vars = n_int
            start_undef = 'INT_MIN'
        elif var_type == 'bool':
            n_vars = n_bool
            start_undef = 'false'

        lines.append('static const cs_' + var_type + '_var_t _' + var_type + '_vars[' + str(n_vars) + '] = {')
        for i, var in enumerate(eval(var_type + '_vars')):
            lines.append('  { // ' + var['name'])
            lines.append('    ' + str(var['index']) + ',')
            lines.append('    ' + str(var['valueReference']) + ',')
            try:
                lines.append('    ' + str(var['causality']) + ',')
            except Exception:
                lines.append('    local,')
            try:
                lines.append('    ' + str(var['variablility']) + ',')
            except Exception:
                lines.append('    continuous,')
            try:
                lines.append('    ' + str(var['initial']) + ',')
            except Exception:
                lines.append('    exact,  // not provided')
            try:
                lines.append('    ' + str(var['start']))
            except Exception:
                lines.append('    ' + start_undef)

            if i+1 < n_vars:
                lines.append('  },')
            else:
                lines.append('  }')
        lines.append('};')
        lines.append('')

    lines.append('static const cs_string_var_t _string_vars[' + str(n_string) + '] = {')
    for i, var in enumerate(string_vars):
        lines.append('  { // ' + var['name'])
        lines.append('    ' + str(var['index']) + ',')
        lines.append('    ' + str(var['valueReference']) + ',')
        try:
            lines.append('    ' + str(var['causality']) + ',')
        except Exception:
            lines.append('    local,')
        try:
            lines.append('    ' + str(var['variablility']) + ',')
        except Exception:
            lines.append('    continuous,')
        try:
            lines.append('    ' + str(var['initial']) + ',')
        except Exception:
            lines.append('    exact,  // not provided')
        try:
            lines.append('    "' + str(var['start']) + '"')
        except Exception:
            lines.append('    ""')

        if i+1 < n_string:
            lines.append('  },')
        else:
            lines.append('  }')
    lines.append('};')
    lines.append('')

    f = open(path, 'w')
    for l in lines:
        f.write(l + os.linesep)
    f.close

#-------------------------------------------------------------------------------

def update_generation_time(path):
    """
    Update generation time in model description
    """

    f = open(path, 'r')

    # Determine or build default names if required

    update = False
    lines = []
    s = 'generationDateAndTime="'
    len_s = len(s)
    for line in f:
        index = line.find(s)
        if index > -1:
            from datetime import datetime
            dt = datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ")
            l = line[0:index] + s + dt + '"'
            line = line[index + len_s:]
            index = line.find('"')
            l += line[index+1:]
            lines.append(l)
            update = True
        else:
            lines.append(line)

    if update:
        f = open(path, 'w')
        for l in lines:
            f.write(l)
        f.close

#-------------------------------------------------------------------------------

if __name__ == '__main__':

    read_from_xml("modelDescription.xml")
    generate_include_file("code_saturne_fmu_variables.h")

    # Update date and time in XML file

    update_generation_time("modelDescription.xml")

    # Next steps needing automation to build FMU from xml:

    # make
    # zip -r code_saturne.fmu modelDescription.xml binaries

    # Note the read_from_xml could be replaces by an alternative method
    # such as reading the tree from the previous fpd, or generating it from a
    # CSV file, in which case the modelDescription.xml file should be
    # generated at the same time, using a a template for the headers
    # and footers.

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
