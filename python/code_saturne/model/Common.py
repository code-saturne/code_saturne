# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2022 EDF S.A.
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
This module defines global constant and the following class and function:
- GuiParam
- dicoLabel
- GuiLabelManager
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import os, logging

#-------------------------------------------------------------------------------
# Global Parameters
#-------------------------------------------------------------------------------

# xml_doc_version modifie le 10/12/07
XML_DOC_VERSION = "2.0"

LABEL_LENGTH_MAX = 32

#-------------------------------------------------------------------------------
# Global GUI parameters
#-------------------------------------------------------------------------------

class GuiParam(object):
    """
    Global options management.
    """
    # 'fr' or 'en' (default)
    #
    try:
        lang = os.environ['LANG'][0:2]
    except Exception:
        lang = 'en'

    # Force English anyway as nearly no translation is available
    lang = 'en'

    # debug
    DEBUG = logging.NOTSET

#-------------------------------------------------------------------------------
# Dictionary : dependance between names and labels
#-------------------------------------------------------------------------------

def dicoLabel(name):
    """
    Correspondence between the names and the labels according to
    whether one is in French or in English.
    """
    for (n, labF, labE) in [('velocity',                      "Vitesse",   "Velocity"),
                            ('pressure',                      "Pression",   "Pressure"),
                            ('hydraulic_head',                "Charge hydraulique", "Hydraulic head"),
                            ('k',                             "EnerTurb",   "TurbEner"),
                            ('epsilon',                       "Dissip",     "Dissip"),
                            ('turbulent_viscosity',           "ViscTurb",   "TurbVisc"),
                            ('r11',                           "R11",        "R11"),
                            ('r22',                           "R22",        "R22"),
                            ('r33',                           "R33",        "R33"),
                            ('r12',                           "R12",        "R12"),
                            ('r13',                           "R13",        "R13"),
                            ('r23',                           "R23",        "R23"),
                            ('rij',                           "Rij",        "Rij"),
                            ('phi',                           "phi",        "phi"),
                            ('alpha',                         "alpha",      "alpha"),
                            ('omega',                         "omega",      "omega"),
                            ('nu_tilda',                      "nu_tilda",   "nu_tilda"),
                            ('smagorinsky_constant^2',        "Csdyn2",     "Csdyn2"),
                            ('temperature_celsius',           "TempC",      "TempC"),
                            ('temperature_kelvin',            "TempK",      "TempK"),
                            ('enthalpy',                      "Enthalpie",  "Enthalpy"),
                            ('potential_temperature',         "TempPot",    "PotTemp"),
                            ('liquid_potential_temperature',  "TempPotLiq", "LiqPotTemp"),
                            ('total_energy',                  "EnerTot",    "TotEner"),
                            ('density',                       "MasseVol",   "Density"),
                            ('molecular_viscosity',           "ViscLam",    "LamVisc"),
                            ('surface_tension',               "TenSurf",    "SurfTen"),
                            ('specific_heat',                 "ChSpec",     "SpecHeat"),
                            ('thermal_conductivity',          "CondTherm",  "ThermalCond"),
                            ('dynamic_diffusion',             "DynDiff",    "DiffDyn"),
                            ('volume_viscosity',              "VolVisc",    "VolVisc"),
                            ('local_time_step',               "pdtlocal",   "LocalTime"),
                            ('courant_number',                "NbCourant",  "CourantNb"),
                            ('fourier_number',                "NbFourier",  "FourierNb"),
                            ('weight_matrix_X',               "VPsolve1",   "VPsolve1"),
                            ('weight_matrix_Y',               "VPsolve2",   "VPsolve2"),
                            ('weight_matrix_Z',               "VPsolve3",   "VPsolve3"),
                            ('est_error_cor_1',               "EsCor1",     "EsCor1"),
                            ('est_error_der_1',               "EsDer1",     "EsDer1"),
                            ('est_error_pre_1',               "EsPre1",     "EsPre1"),
                            ('est_error_tot_1',               "EsTot1",     "EsTot1"),
                            ('est_error_cor_2',               "EsCor2",     "EsCor2"),
                            ('est_error_der_2',               "EsDer2",     "EsDer2"),
                            ('est_error_pre_2',               "EsPre2",     "EsPre2"),
                            ('est_error_tot_2',               "EsTot2",     "EsTot2")]:

        if n == name:
            if GuiParam.lang == 'fr':
                label = labF
            else:
                label = labE

    return label

#-------------------------------------------------------------------------------
# Gui Parameters for autocompletion and forbidden labels
#-------------------------------------------------------------------------------

class GuiLabelManager(object):

    def __init__(self):

        self._completer = {}
        self._completer['mesh_selection'] = ['all[]',
                                            'plane[a, b, c, d, epsilon]',
                                            'plane[a, b, c, d, inside]',
                                            'plane[a, b, c, d, outside]',
                                            'plane[n_x, n_y, n_z, x0, y0, z0, epsilon]',
                                            'plane[n_x, n_y, n_z, x0, y0, z0, inside]',
                                            'plane[n_x, n_y, n_z, x0, y0, z0, outside]',
                                            'box[xmin, ymin, zmin, xmax, ymax, zmax]',
                                            'box[x0, y0, z0, dx1, dy1, dz1, dx2, dy2, dz2, dx3, dy3, dz3]',
                                            'cylinder[x0, y0, z0, x1, y1, z1, radius]',
                                            'sphere[x_c, y_c, z_c, radius]']

        self._completer['mesh_normal'] = ['normal[x, y, z, epsilon]']

        self._forbidden = {}
        self._forbidden['notebook'] = ['temperature',
                                       'pressure',
                                       'density',
                                       'specific_heat',
                                       'molecular_viscosty',
                                       'thermal_conductivity',
                                       'rho', 'rho0', 'ro0',
                                       'mu', 'mu0', 'viscl0',
                                       'p0', 'cp0', 'cp', 'lambda0',
                                       'enthalpy', 'volume_fraction',
                                       'vol_f', 'uref', 'almax']


    def getCompleter(self, key):

        if key not in self._completer.keys():
            raise Exception("Uknown key: %s" % key)

        return self._completer[key]

    def getForbidden(self, key):

        if key not in self._forbidden.keys():
            raise Exception("Uknown key: %s" % key)

        return self._forbidden[key]

#-------------------------------------------------------------------------------
# End of Common
#-------------------------------------------------------------------------------
