# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2020 EDF S.A.
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
# End of Common
#-------------------------------------------------------------------------------
