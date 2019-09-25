# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2019 EDF S.A.
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
This module defines the following function:
- displaySelectedPage
"""

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.model.Common import *

#-------------------------------------------------------------------------------
# displaySelectedPage direct to the good page with its name
#-------------------------------------------------------------------------------

def displaySelectedPage(page_name, root, case, stbar=None, tree=None):
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
        thisPage = Page.IdentityAndPathesView(root, case)

    elif page_name == tr("Mesh"):
        import code_saturne.Pages.SolutionDomainView as Page
        thisPage = Page.SolutionDomainView(root, case, stbar, tree)

    elif page_name == tr("Preprocessing"):
        import code_saturne.Pages.PreprocessingView as Page
        thisPage = Page.PreprocessingView(root, case, stbar)

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
        import code_saturne.Pages.ThermalView as Page
        thisPage = Page.ThermalView(root, case, tree)

    elif page_name == tr("Gas combustion"):
        import code_saturne.Pages.GasCombustionView as Page
        thisPage = Page.GasCombustionView(root, case)

    elif page_name == tr("Pulverized fuel combustion"):
        import code_saturne.Pages.CoalCombustionView as Page
        thisPage = Page.CoalCombustionView(root, case, stbar)

    elif page_name == tr("Electrical models"):
        import code_saturne.Pages.ElectricalView as Page
        thisPage = Page.ElectricalView(root, case, stbar)

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
        if case.xmlRootNode().tagName == "NEPTUNE_CFD_GUI":
            import code_saturne.Pages.MainFieldsSourceTermsView as Page
            thisPage = Page.MainFieldsSourceTermsView(root, case, stbar)
        else:
            import code_saturne.Pages.SourceTermsView as Page
            thisPage = Page.SourceTermsView(root, case, stbar)

    elif page_name == tr("Groundwater laws"):
        import code_saturne.Pages.GroundwaterLawView as Page
        thisPage = Page.GroundwaterLawView(root, case)

    elif page_name == tr("Fluid properties"):
        import code_saturne.Pages.FluidCharacteristicsView as Page
        thisPage = Page.FluidCharacteristicsView(root, case)

    elif page_name == tr("Body forces"):
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

    elif page_name == tr("Start/Restart"):
        import code_saturne.Pages.StartRestartView as Page
        thisPage = Page.StartRestartView(root, case)

    elif page_name == tr("Postprocessing"):
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

    elif page_name == tr("Performance settings"):
        import code_saturne.Pages.PerformanceTuningView as Page
        thisPage = Page.PerformanceTuningView(root, case)

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
# End of Toolbox
#-------------------------------------------------------------------------------
