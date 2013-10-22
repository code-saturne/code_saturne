# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2013 EDF S.A.
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
This module defines the following classes and functions:
- GuiParam
- displaySelectedPage
- dicoLabel
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import os, sys, logging

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from Base.Common import *

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
    #
    DEBUG = logging.NOTSET

#-------------------------------------------------------------------------------
# displaySelectedPage direct to the good page with its name
#-------------------------------------------------------------------------------

def displaySelectedPage(page_name, root, case, stbar=None, study=None, tree=None):
    """
    This function enables to display a new page when the TreeNavigator
    send the order.
    """
    # 'win' is the frame-support of the Pages
    # 'thisPage' is the instance of classes which create thePages
    # 'page_name' is the name of the page
    #

    if page_name == tr("Identity and paths"):
        import Pages.IdentityAndPathesView as Page
        thisPage = Page.IdentityAndPathesView(root, case, study)

    elif page_name == tr("Meshes selection"):
        import Pages.SolutionDomainView as Page
        thisPage = Page.SolutionDomainView(root, case, stbar)

    elif page_name == tr("Mesh quality criteria"):
        import Pages.SolutionVerifView as Page
        thisPage = Page.SolutionVerifView(root, case)

    elif page_name == tr("Volume regions definition"):
        import Pages.LocalizationView as Page
        thisPage = Page.VolumeLocalizationView(root, case, tree)

    elif page_name == tr("Calculation features"):
        import Pages.AnalysisFeaturesView as Page
        thisPage = Page.AnalysisFeaturesView(root, case, tree)

    elif page_name == tr("Deformable mesh"):
        import Pages.MobileMeshView as Page
        thisPage = Page.MobileMeshView(root, case, tree)

    elif page_name == tr("Turbulence models"):
        import Pages.TurbulenceView as Page
        thisPage = Page.TurbulenceView(root, case)

    elif page_name == tr("Thermal model"):
        import Pages.ThermalScalarView as Page
        thisPage = Page.ThermalScalarView(root, case, tree)

    elif page_name == tr("Gas combustion"):
        import Pages.GasCombustionView as Page
        thisPage = Page.GasCombustionView(root, case)

    elif page_name == tr("Pulverized fuel combustion"):
        import Pages.CoalCombustionView as Page
        thisPage = Page.CoalCombustionView(root, case, stbar)

    elif page_name == tr("Electrical models"):
        import Pages.ElectricalView as Page
        thisPage = Page.ElectricalView(root, case, stbar)

    elif page_name == tr("Radiative transfers"):
        import Pages.ThermalRadiationView as Page
        thisPage = Page.ThermalRadiationView(root, case, tree)

    elif page_name == tr("Conjugate heat transfer"):
        import Pages.ConjugateHeatTransferView as Page
        thisPage = Page.ConjugateHeatTransferView(root, case)

    elif page_name == tr("Initialization"):
        import Pages.InitializationView as Page
        thisPage = Page.InitializationView(root, case, stbar)

    elif page_name == tr("Head losses"):
        import Pages.HeadLossesView as Page
        thisPage = Page.HeadLossesView(root, case)

    elif page_name == tr("Source terms"):
        import Pages.SourceTermsView as Page
        thisPage = Page.SourceTermsView(root, case, stbar)

    elif page_name == tr("Coriolis Source Terms"):
        import Pages.CoriolisSourceTermsView as Page
        thisPage = Page.CoriolisSourceTermsView(root, case)

    elif page_name == tr("Reference values"):
        import Pages.ReferenceValuesView as Page
        thisPage = Page.ReferenceValuesView(root, case)

    elif page_name == tr("Fluid properties"):
        import Pages.FluidCharacteristicsView as Page
        thisPage = Page.FluidCharacteristicsView(root, case)

    elif page_name == tr("Gravity"):
        import Pages.BodyForcesView as Page
        thisPage = Page.BodyForcesView(root, case)

    elif page_name == tr("Species transport"):
        import Pages.DefineUserScalarsView as Page
        thisPage = Page.DefineUserScalarsView(root, case, stbar)

    elif page_name == tr("Global settings"):
        import Pages.LagrangianView as Page
        thisPage = Page.LagrangianView(root, case)

    elif page_name == tr("Statistics"):
        import Pages.LagrangianStatisticsView as Page
        thisPage = Page.LagrangianStatisticsView(root, case)

    elif page_name == tr("Definition of boundary regions"):
        import Pages.LocalizationView as Page
        thisPage = Page.BoundaryLocalizationView(root, case, tree)

    elif page_name == tr("Boundary conditions"):
        import Pages.BoundaryConditionsView as Page
        thisPage = Page.BoundaryConditionsView(root, case)

    elif page_name == tr("Particles boundary conditions"):
        import Pages.LagrangianBoundariesView as Page
        thisPage = Page.LagrangianBoundariesView(root, case)

    elif page_name == tr("Time averages"):
        import Pages.TimeAveragesView as Page
        thisPage = Page.TimeAveragesView(root, case, stbar)

    elif page_name == tr("Time step"):
        import Pages.TimeStepView as Page
        thisPage = Page.TimeStepView(root, case)

    elif page_name == tr("Pseudo-Time step"):
        import Pages.TimeStepView as Page
        thisPage = Page.TimeStepView(root, case)

    elif page_name == tr("Steady flow management"):
        import Pages.SteadyManagementView as Page
        thisPage = Page.SteadyManagementView(root, case)

    elif page_name == tr("Output control"):
        import Pages.OutputControlView as Page
        thisPage = Page.OutputControlView(root, case, tree)

    elif page_name == tr("Volume solution control"):
        import Pages.OutputVolumicVariablesView as Page
        thisPage = Page.OutputVolumicVariablesView(root, case)

    elif page_name == tr("Surface solution control"):
        import Pages.OutputSurfacicVariablesView as Page
        thisPage = Page.OutputSurfacicVariablesView(root, case)

    elif page_name == tr("Lagrangian solution control"):
        import Pages.LagrangianOutputView as Page
        thisPage = Page.LagrangianOutputView(root, case)

    elif page_name == tr("Profiles"):
        import Pages.ProfilesView as Page
        thisPage = Page.ProfilesView(root, case, stbar)

    elif page_name == tr("Equation parameters"):
        import Pages.NumericalParamEquationView as Page
        thisPage = Page.NumericalParamEquationView(root, case)

    elif page_name == tr("Global parameters"):
        import Pages.NumericalParamGlobalView as Page
        thisPage = Page.NumericalParamGlobalView(root, case, tree)

    elif page_name == tr("Start/Restart"):
        import Pages.StartRestartView as Page
        thisPage = Page.StartRestartView(root, case)

    elif page_name == tr("Performance tuning"):
        import Pages.PerformanceTuningView as Page
        thisPage = Page.PerformanceTuningView(root, case)

    elif page_name == tr("Prepare batch calculation"):
        import Pages.BatchRunningView as Page
        thisPage = Page.BatchRunningView(root, case)

    elif page_name == tr("Fluid structure interaction"):
        import Pages.FluidStructureInteractionView as Page
        thisPage = Page.FluidStructureInteractionView(root, case)

    elif page_name == tr("Atmospheric flows"):
        import Pages.AtmosphericFlowsView as Page
        thisPage = Page.AtmosphericFlowsView(root, case)


    else:
        msg = tr("Warning: the corresponding Page %s doesn't exist!") % page_name
        print(msg)
        # So we display the Welcome Page!
        import Pages.WelcomeView as Page
        thisPage = Page.WelcomeView()

    case['current_page'] = str(page_name)

    return thisPage


def tr(text):
    """
    Translation
    """
    return text

#-------------------------------------------------------------------------------
# Dictionary : dependance between names and labels
#-------------------------------------------------------------------------------

def dicoLabel(name):
    """
    Correspondence between the names and the labels according to
    whether one is in French or in English.
    """
    for (n, labF, labE) in [('velocity_U',                    "VitesseX",   "VelocityX"),
                            ('velocity_V',                    "VitesseY",   "VelocityY"),
                            ('velocity_W',                    "VitesseZ",   "VelocityZ"),
                            ('pressure',                      "Pression",   "Pressure"),
                            ('turb_k',                        "EnerTurb",   "TurbEner"),
                            ('turb_eps',                      "Dissip",     "Dissip"),
                            ('turb_viscosity',                "ViscTurb",   "TurbVisc"),
                            ('component_R11',                 "R11",        "R11"),
                            ('component_R22',                 "R22",        "R22"),
                            ('component_R33',                 "R33",        "R33"),
                            ('component_R12',                 "R12",        "R12"),
                            ('component_R13',                 "R13",        "R13"),
                            ('component_R23',                 "R23",        "R23"),
                            ('turb_phi',                      "phi",        "phi"),
                            ('turb_alpha',                    "alpha",      "alpha"),
                            ('turb_omega',                    "omega",      "omega"),
                            ('turb_nusa',                     "nusa",       "nusa"),
                            ('smagorinsky_constant',          "Csdyn2",     "Csdyn2"),
                            ('temperature_celsius',           "TempC",      "TempC"),
                            ('temperature_kelvin',            "TempK",      "TempK"),
                            ('enthalpy',                      "Enthalpie",  "Enthalpy"),
                            ('potential_temperature',         "TempPot",    "PotTemp"),
                            ('liquid_potential_temperature',  "TempPotLiq", "LiqPotTemp"),
                            ('density',                       "MasseVol",   "Density"),
                            ('molecular_viscosity',           "ViscLam",    "LamVisc"),
                            ('specific_heat',                 "ChSpec",     "SpecHeat"),
                            ('thermal_conductivity',          "CondTherm",  "ThermalCond"),
                            ('dynamic_diffusion',             "DynDiff",    "DiffDyn"),
                            ('volumic_viscosity',             "VolVisc",    "VolVisc"),
                            ('local_time_step',               "pdtlocal",   "LocalTime"),
                            ('courant_number',                "NbCourant",  "CourantNb"),
                            ('fourier_number',                "NbFourier",  "FourierNb"),
                            ('weight_matrix_X',               "VPsolve1",   "VPsolve1"),
                            ('weight_matrix_Y',               "VPsolve2",   "VPsolve2"),
                            ('weight_matrix_Z',               "VPsolve3",   "VPsolve3")]:

        if n == name:
            if GuiParam.lang == 'fr':
                label = labF
            else:
                label = labE

    return label

#-------------------------------------------------------------------------------
# End of Toolbox
#-------------------------------------------------------------------------------
