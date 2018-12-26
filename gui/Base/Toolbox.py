# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2018 EDF S.A.
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

from code_saturne.Base.Common import *

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

    if page_name == tr("Calculation environment"):
        import code_saturne.Pages.IdentityAndPathesView as Page
        thisPage = Page.IdentityAndPathesView(root, case, study)

    elif page_name == tr("Meshes selection"):
        import code_saturne.Pages.SolutionDomainView as Page
        thisPage = Page.SolutionDomainView(root, case, stbar)

    elif page_name == tr("Notebook"):
        import code_saturne.Pages.NotebookView as Page
        thisPage = Page.NotebookView(root, case)

    elif page_name == tr("Volume zones"):
        import code_saturne.Pages.LocalizationView as Page
        thisPage = Page.VolumeLocalizationView(root, case, tree)

    elif page_name == tr("Calculation features"):
        import code_saturne.Pages.AnalysisFeaturesView as Page
        thisPage = Page.AnalysisFeaturesView(root, case, tree)

    elif page_name == tr("Deformable mesh"):
        import code_saturne.Pages.MobileMeshView as Page
        thisPage = Page.MobileMeshView(root, case, tree)

    elif page_name == tr("Turbulence models"):
         if case.xmlRootNode().tagName == "NEPTUNE_CFD_GUI" :
            import code_saturne.Pages.TurbulenceNeptuneView as Page
            thisPage = Page.TurbulenceView(root, case)
         else :
            import code_saturne.Pages.TurbulenceView as Page
            thisPage = Page.TurbulenceView(root, case)

    elif page_name == tr("Thermal model"):
        import code_saturne.Pages.ThermalScalarView as Page
        thisPage = Page.ThermalScalarView(root, case, tree)

    elif page_name == tr("Gas combustion"):
        import code_saturne.Pages.GasCombustionView as Page
        thisPage = Page.GasCombustionView(root, case)

    elif page_name == tr("Pulverized fuel combustion"):
        import code_saturne.Pages.CoalCombustionView as Page
        thisPage = Page.CoalCombustionView(root, case, stbar)

    elif page_name == tr("Electrical models"):
        import code_saturne.Pages.ElectricalView as Page
        thisPage = Page.ElectricalView(root, case, stbar)

    elif page_name == tr("Radiative transfers"):
        import code_saturne.Pages.ThermalRadiationView as Page
        thisPage = Page.ThermalRadiationView(root, case, tree)

    elif page_name == tr("Conjugate heat transfer"):
        import code_saturne.Pages.ConjugateHeatTransferView as Page
        thisPage = Page.ConjugateHeatTransferView(root, case)

    elif page_name == tr("Main fields initialization"):
        import code_saturne.Pages.MainFieldsInitializationView as Page
        thisPage = Page.MainFieldsInitializationView(root, case)

    elif page_name == tr("Initialization"):
        import code_saturne.Pages.InitializationView as Page
        thisPage = Page.InitializationView(root, case, stbar)

    elif page_name == tr("Head losses"):
        import code_saturne.Pages.HeadLossesView as Page
        thisPage = Page.HeadLossesView(root, case)

    elif page_name == tr("Porosity"):
        import code_saturne.Pages.PorosityView as Page
        thisPage = Page.PorosityView(root, case)

    elif page_name == tr("Source terms"):
        import code_saturne.Pages.SourceTermsView as Page
        thisPage = Page.SourceTermsView(root, case, stbar)

    elif page_name == tr("Coriolis Source Terms"):
        import code_saturne.Pages.CoriolisSourceTermsView as Page
        thisPage = Page.CoriolisSourceTermsView(root, case)

    elif page_name == tr("Groundwater laws"):
        import code_saturne.Pages.GroundwaterLawView as Page
        thisPage = Page.GroundwaterLawView(root, case)

    elif page_name == tr("Physical properties"):
        import code_saturne.Pages.ReferenceValuesView as Page
        thisPage = Page.ReferenceValuesView(root, case)

    elif page_name == tr("Fluid properties"):
        import code_saturne.Pages.FluidCharacteristicsView as Page
        thisPage = Page.FluidCharacteristicsView(root, case)

    elif page_name == tr("Gravity"):
        import code_saturne.Pages.BodyForcesView as Page
        thisPage = Page.BodyForcesView(root, case)

    elif page_name == tr("Species transport"):
        if case.xmlRootNode().tagName == "NEPTUNE_CFD_GUI" :
            import code_saturne.Pages.SpeciesView as Page
            thisPage = Page.SpeciesView(root, case)
        else :
            import code_saturne.Pages.DefineUserScalarsView as Page
            thisPage = Page.DefineUserScalarsView(root, case, stbar, tree)

    elif page_name == tr("Turbomachinery"):
        import code_saturne.Pages.TurboMachineryView as Page
        thisPage = Page.TurboMachineryView(root, case)

    elif page_name == tr("Fans"):
        import code_saturne.Pages.FansView as Page
        thisPage = Page.FansView(root, case)

    elif page_name == tr("Groundwater flows"):
        import code_saturne.Pages.GroundwaterView as Page
        thisPage = Page.GroundwaterView(root, case)

    elif page_name == tr("Particles and droplets tracking"):
        import code_saturne.Pages.LagrangianView as Page
        thisPage = Page.LagrangianView(root, case)

    elif page_name == tr("Statistics"):
        import code_saturne.Pages.LagrangianStatisticsView as Page
        thisPage = Page.LagrangianStatisticsView(root, case)

    elif page_name == tr("Main fields boundary conditions"):
        import code_saturne.Pages.BoundaryConditionsViewNeptune as Page
        thisPage = Page.BoundaryConditionsView(root, case)

    elif page_name == tr("Cathare Coupling"):
        import code_saturne.Pages.CathareCouplingView as Page
        thisPage = Page.CathareCouplingView(root, case)

    elif page_name == tr("Boundary zones"):
        import code_saturne.Pages.LocalizationView as Page
        thisPage = Page.BoundaryLocalizationView(root, case, tree)

    elif page_name == tr("Boundary conditions"):
        if case.xmlRootNode().tagName == "NEPTUNE_CFD_GUI":
            import code_saturne.Pages.BoundaryConditionsViewNeptune as Page
            thisPage = Page.BoundaryConditionsView(root, case)
        else:
            import code_saturne.Pages.BoundaryConditionsView as Page
            thisPage = Page.BoundaryConditionsView(root, case)

    elif page_name == tr("Particle boundary conditions"):
        import code_saturne.Pages.LagrangianBoundariesView as Page
        thisPage = Page.LagrangianBoundariesView(root, case)

    elif page_name == tr("Time averages"):
        import code_saturne.Pages.TimeAveragesView as Page
        thisPage = Page.TimeAveragesView(root, case, stbar)

    elif page_name == tr("Time settings"):
        if case.xmlRootNode().tagName == "NEPTUNE_CFD_GUI" :
            import code_saturne.Pages.TimeStepViewNeptune as Page
            thisPage = Page.TimeStepView(root, case)
        else:
            import code_saturne.Pages.TimeStepView as Page
            thisPage = Page.TimeStepView(root, case, tree)

    elif page_name == tr("Postprocessing"):
        if case.xmlRootNode().tagName == "NEPTUNE_CFD_GUI" :
            import code_saturne.Pages.OutputControlViewNeptune as Page
            thisPage = Page.OutputControlView(root, case, tree)
        else:
            import code_saturne.Pages.OutputControlView as Page
            thisPage = Page.OutputControlView(root, case, tree)

    elif page_name == tr("Additional user arrays"):
        import code_saturne.Pages.UsersControlView as Page
        thisPage = Page.UsersControlView(root, case)

    elif page_name == tr("Volume solution control"):
        import code_saturne.Pages.OutputVolumicVariablesView as Page
        thisPage = Page.OutputVolumicVariablesView(root, case)

    elif page_name == tr("Surface solution control"):
        import code_saturne.Pages.OutputSurfacicVariablesView as Page
        thisPage = Page.OutputSurfacicVariablesView(root, case)

    elif page_name == tr("Lagrangian solution control"):
        import code_saturne.Pages.LagrangianOutputView as Page
        thisPage = Page.LagrangianOutputView(root, case)

    elif page_name == tr("Profiles"):
        import code_saturne.Pages.ProfilesView as Page
        thisPage = Page.ProfilesView(root, case, stbar)

    elif page_name == tr("Balance by zone"):
        import code_saturne.Pages.BalanceView as Page
        thisPage = Page.BalanceView(root, case)

    elif page_name == tr("Equation parameters"):
        if case.xmlRootNode().tagName == "NEPTUNE_CFD_GUI" :
            import code_saturne.Pages.NumericalParamEquationViewNeptune as Page
            thisPage = Page.NumericalParamEquationView(root, case)
        else :
            import code_saturne.Pages.NumericalParamEquationView as Page
            thisPage = Page.NumericalParamEquationView(root, case)

    elif page_name == tr("Numerical parameters"):
        if case.xmlRootNode().tagName == "NEPTUNE_CFD_GUI" :
            import code_saturne.Pages.GlobalNumericalParametersView as Page
            thisPage = Page.GlobalNumericalParametersView(root, case)
        else :
            import code_saturne.Pages.NumericalParamGlobalView as Page
            thisPage = Page.NumericalParamGlobalView(root, case, tree)

    elif page_name == tr("Calculation management"):
        import code_saturne.Pages.StartRestartView as Page
        thisPage = Page.StartRestartView(root, case)

    elif page_name == tr("Performance tuning"):
        import code_saturne.Pages.PerformanceTuningView as Page
        thisPage = Page.PerformanceTuningView(root, case)

    elif page_name == tr("Prepare batch calculation"):
        import code_saturne.Pages.BatchRunningView as Page
        thisPage = Page.BatchRunningView(root, case)

    elif page_name == tr("Fluid structure interaction"):
        import code_saturne.Pages.FluidStructureInteractionView as Page
        thisPage = Page.FluidStructureInteractionView(root, case)

    elif page_name == tr("Atmospheric flows"):
        import code_saturne.Pages.AtmosphericFlowsView as Page
        thisPage = Page.AtmosphericFlowsView(root, case)

    elif page_name == tr("Non condensable gases"):
        import code_saturne.Pages.NonCondensableView as Page
        thisPage = Page.NonCondensableView(root, case, tree)

    elif page_name == tr("Main fields"):
        import code_saturne.Pages.MainFieldsView as Page
        thisPage = Page.MainFieldsView(root, case, tree)

    elif page_name == tr("Thermodynamics"):
        import code_saturne.Pages.ThermodynamicsView as Page
        thisPage = Page.ThermodynamicsView(root, case)

    elif page_name == tr("Closure modeling"):
        import code_saturne.Pages.InterfacialForcesView as Page
        thisPage = Page.InterfacialForcesView(root, case, tree)

    elif page_name == tr("Interfacial enthalpy transfer"):
        import code_saturne.Pages.InterfacialEnthalpyView as Page
        thisPage = Page.InterfacialEnthalpyView(root, case)

    elif page_name == tr("Nucleate boiling parameters"):
        import code_saturne.Pages.NucleateBoilingView as Page
        thisPage = Page.NucleateBoilingView(root, case)

    elif page_name == tr("Droplet condensation-evaporation"):
        import code_saturne.Pages.DropletCondensationEvaporationView as Page
        thisPage = Page.DropletCondensationEvaporationView(root, case)

    elif page_name == tr("Particles interactions"):
        import code_saturne.Pages.SolidView as Page
        thisPage = Page.SolidView(root, case)

    elif page_name == tr("Interfacial area"):
        import code_saturne.Pages.InterfacialAreaView as Page
        thisPage = Page.InterfacialAreaView(root, case)

    elif page_name == tr("OpenTurns study"):
        import code_saturne.Pages.OpenTurnsView as Page
        thisPage = Page.OpenTurnsView(root, case)

    else:
        msg = tr("Warning: the corresponding Page %s doesn't exist!") % page_name
        print(msg)
        # So we display the Welcome Page!
        import code_saturne.Pages.WelcomeView as Page
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
# End of Toolbox
#-------------------------------------------------------------------------------
