# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2023 EDF S.A.
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

from code_saturne.gui.base.QtCore import QCoreApplication
from code_saturne.model.Common import *


#-------------------------------------------------------------------------------
# displaySelectedPage direct to the good page with its name
#-------------------------------------------------------------------------------

class NonExistentPage(Exception):

    def __init__(self, page_name):
        msg = tr("Warning: the corresponding Page %s doesn't exist!") % page_name
        super().__init__(msg)


def displaySelectedPage(page_name, root, case, stbar=None, tree=None):
    """
    This function enables to display a new page when the TreeNavigator
    send the order.
    """
    # 'win' is the frame-support of the Pages
    # 'thisPage' is the instance of classes which create thePages
    # 'page_name' is the name of the page
    #

    try:
        thisPage = displayStaticPage(case, page_name, root, stbar, tree)
    except NonExistentPage:
        index = case['current_index']
        item = index.internalPointer()
        zone_name = item.itemData[0]

        if item.parentItem.itemData[0] == tr("Boundary conditions"):
            if case.xmlRootNode().tagName == tr("NEPTUNE_CFD_GUI"):
                import code_saturne.gui.case.BoundaryConditionsViewNeptune as Page
                thisPage = Page.BoundaryConditionsView(root, case, zone_name)
            else:
                import code_saturne.gui.case.BoundaryConditionsView as Page
                thisPage = Page.BoundaryConditionsView(root, case, zone_name)

        elif item.parentItem.itemData[0] == tr("Volume conditions"):
            import code_saturne.gui.case.VolumicConditionsView as Page
            thisPage = Page.VolumicConditionsView(root, case, zone_name)

        elif item.parentItem.itemData[0] == tr("Volume zones"):
            import code_saturne.gui.case.LocalizationView as Page
            thisPage = Page.VolumeLocalizationView(root, case, tree, hide_all=True)

        else:
            import code_saturne.gui.case.WelcomeView as Page
            thisPage = Page.WelcomeView()

    case['current_page'] = str(page_name)

    return thisPage


def displayStaticPage(case, page_name, root, stbar, tree):
    if page_name == tr("Calculation environment"):
        import code_saturne.gui.case.IdentityAndPathesView as Page
        thisPage = Page.IdentityAndPathesView(root, case)

    elif page_name == tr("Mesh"):
        import code_saturne.gui.case.SolutionDomainView as Page
        thisPage = Page.SolutionDomainView(root, case, stbar, tree)

    elif page_name == tr("Preprocessing"):
        import code_saturne.gui.case.PreprocessingView as Page
        thisPage = Page.PreprocessingView(root, case, stbar)

    elif page_name == tr("Volume zones"):
        import code_saturne.gui.case.LocalizationView as Page
        thisPage = Page.VolumeLocalizationView(root, case, tree)

    elif page_name == tr("Boundary zones"):
        import code_saturne.gui.case.LocalizationView as Page
        thisPage = Page.BoundaryLocalizationView(root, case, tree)

    elif page_name == tr("Notebook"):
        import code_saturne.gui.case.NotebookView as Page
        thisPage = Page.NotebookView(root, case)

    elif page_name == tr("Time tables"):
        import code_saturne.gui.case.TimeTablesView as Page
        thisPage = Page.TimeTablesView(root, case, stbar, tree)

    elif page_name == tr("Volume conditions"):
        import code_saturne.gui.case.VolumicNatureView as Page
        thisPage = Page.VolumicNatureView(root, case, tree)

    elif page_name == tr("Calculation features"):
        import code_saturne.gui.case.AnalysisFeaturesView as Page
        thisPage = Page.AnalysisFeaturesView(root, case, tree)

    elif page_name == tr("Deformable mesh"):
        import code_saturne.gui.case.MobileMeshView as Page
        thisPage = Page.MobileMeshView(root, case, tree)

    elif page_name == tr("Gas combustion"):
        import code_saturne.gui.case.GasCombustionView as Page
        thisPage = Page.GasCombustionView(root, case)

    elif page_name == tr("Pulverized fuel combustion"):
        import code_saturne.gui.case.CoalCombustionView as Page
        thisPage = Page.CoalCombustionView(root, case, stbar)

    elif page_name == tr("Electrical models"):
        import code_saturne.gui.case.ElectricalView as Page
        thisPage = Page.ElectricalView(root, case, stbar)

    elif page_name == tr("Groundwater flows"):
        import code_saturne.gui.case.GroundwaterView as Page
        thisPage = Page.GroundwaterView(root, case)

    elif page_name == tr("Thermal model"):
        import code_saturne.gui.case.ThermalView as Page
        thisPage = Page.ThermalView(root, case, tree)

    elif page_name == tr("Turbulence models"):
        if case.xmlRootNode().tagName == "NEPTUNE_CFD_GUI":
            import code_saturne.gui.case.TurbulenceNeptuneView as Page
            thisPage = Page.TurbulenceView(root, case)
        else:
            import code_saturne.gui.case.TurbulenceView as Page
            thisPage = Page.TurbulenceView(root, case)

    elif page_name == tr("Body forces"):
        import code_saturne.gui.case.BodyForcesView as Page
        thisPage = Page.BodyForcesView(root, case)

    elif page_name == tr("Species transport"):
        if case.xmlRootNode().tagName == "NEPTUNE_CFD_GUI":
            import code_saturne.gui.case.SpeciesView as Page
            thisPage = Page.SpeciesView(root, case)
        else:
            import code_saturne.gui.case.DefineUserScalarsView as Page
            thisPage = Page.DefineUserScalarsView(root, case, stbar, tree)

    elif page_name == tr("Turbomachinery"):
        import code_saturne.gui.case.TurboMachineryView as Page
        thisPage = Page.TurboMachineryView(root, case)

    elif page_name == tr("Fans"):
        import code_saturne.gui.case.FansView as Page
        thisPage = Page.FansView(root, case)

    elif page_name == tr("Particles and droplets tracking"):
        import code_saturne.gui.case.LagrangianView as Page
        thisPage = Page.LagrangianView(root, case)

    elif page_name == tr("Statistics"):
        import code_saturne.gui.case.LagrangianStatisticsView as Page
        thisPage = Page.LagrangianStatisticsView(root, case)

    elif page_name == tr("Main fields boundary conditions"):
        import code_saturne.gui.case.BoundaryConditionsViewNeptune as Page
        thisPage = Page.BoundaryConditionsView(root, case)

    elif page_name == tr("Immersed Boundaries"):
        import code_saturne.gui.case.ImmersedBoundariesViewNeptune as Page
        thisPage = Page.ImmersedBoundariesViewNeptune(root, case)

    elif page_name == tr("Boundary conditions"):
        import code_saturne.gui.case.BoundaryNatureView as Page
        thisPage = Page.BoundaryNatureView(root, case, tree)

    elif page_name == tr("Coupling parameters"):
        import code_saturne.gui.case.CouplingParametersView as Page
        thisPage = Page.CouplingParametersView(root, case)

    elif page_name == tr("Time averages"):
        import code_saturne.gui.case.TimeAveragesView as Page
        thisPage = Page.TimeAveragesView(root, case, stbar)

    elif page_name == tr("Time settings"):
        if case.xmlRootNode().tagName == "NEPTUNE_CFD_GUI":
            import code_saturne.gui.case.TimeStepViewNeptune as Page
            thisPage = Page.TimeStepView(root, case)
        else:
            import code_saturne.gui.case.TimeStepView as Page
            thisPage = Page.TimeStepView(root, case, tree)

    elif page_name == tr("Start/Restart"):
        import code_saturne.gui.case.StartRestartView as Page
        thisPage = Page.StartRestartView(root, case)

    elif page_name == tr("Postprocessing"):
        import code_saturne.gui.case.OutputControlView as Page
        thisPage = Page.OutputControlView(root, case, tree)

    elif page_name == tr("Calculator"):
        import code_saturne.gui.case.UserCalculatorView as Page
        thisPage = Page.UserCalculatorView(root, case)

    elif page_name == tr("Additional user arrays"):
        import code_saturne.gui.case.UsersControlView as Page
        thisPage = Page.UsersControlView(root, case)

    elif page_name == tr("Volume solution control"):
        import code_saturne.gui.case.OutputVolumicVariablesView as Page
        thisPage = Page.OutputVolumicVariablesView(root, case)

    elif page_name == tr("Surface solution control"):
        import code_saturne.gui.case.OutputSurfacicVariablesView as Page
        thisPage = Page.OutputSurfacicVariablesView(root, case)

    elif page_name == tr("Lagrangian solution control"):
        import code_saturne.gui.case.LagrangianOutputView as Page
        thisPage = Page.LagrangianOutputView(root, case)

    elif page_name == tr("Profiles"):
        import code_saturne.gui.case.ProfilesView as Page
        thisPage = Page.ProfilesView(root, case, stbar)

    elif page_name == tr("Balance by zone"):
        import code_saturne.gui.case.BalanceView as Page
        thisPage = Page.BalanceView(root, case)

    elif page_name == tr("Equation parameters"):
        if case.xmlRootNode().tagName == "NEPTUNE_CFD_GUI":
            import code_saturne.gui.case.NumericalParamEquationViewNeptune as Page
            thisPage = Page.NumericalParamEquationView(root, case)
        else:
            import code_saturne.gui.case.NumericalParamEquationView as Page
            thisPage = Page.NumericalParamEquationView(root, case)

    elif page_name == tr("Numerical parameters"):
        if case.xmlRootNode().tagName == "NEPTUNE_CFD_GUI":
            import code_saturne.gui.case.GlobalNumericalParametersView as Page
            thisPage = Page.GlobalNumericalParametersView(root, case)
        else:
            import code_saturne.gui.case.NumericalParamGlobalView as Page
            thisPage = Page.NumericalParamGlobalView(root, case, tree)

    elif page_name == tr("Performance settings"):
        import code_saturne.gui.case.PerformanceTuningView as Page
        thisPage = Page.PerformanceTuningView(root, case)

    elif page_name == tr("Atmospheric flows"):
        import code_saturne.gui.case.AtmosphericFlowsView as Page
        thisPage = Page.AtmosphericFlowsView(root, case)

    elif page_name == tr("Non condensable gases"):
        import code_saturne.gui.case.NonCondensableView as Page
        thisPage = Page.NonCondensableView(root, case, tree)

    elif page_name == tr("Main fields"):
        import code_saturne.gui.case.MainFieldsView as Page
        thisPage = Page.MainFieldsView(root, case, tree)

    elif page_name == tr("Closure modeling"):
        import code_saturne.gui.case.InterfacialForcesView as Page
        thisPage = Page.InterfacialForcesView(root, case, tree)

    elif page_name == tr("Interfacial enthalpy transfer"):
        import code_saturne.gui.case.InterfacialEnthalpyView as Page
        thisPage = Page.InterfacialEnthalpyView(root, case)

    elif page_name == tr("Wall transfer parameters"):
        import code_saturne.gui.case.NeptuneWallTransferView as Page
        thisPage = Page.NeptuneWallTransferView(root, case)

    elif page_name == tr("Particles interactions"):
        import code_saturne.gui.case.SolidView as Page
        thisPage = Page.SolidView(root, case)

    elif page_name == tr("Interfacial area"):
        import code_saturne.gui.case.InterfacialAreaView as Page
        thisPage = Page.InterfacialAreaView(root, case)

    else:
        raise NonExistentPage(page_name)
        ## So we display the Welcome Page!


    return thisPage


def tr(text):
    """
    Translation
    """

    # Note that the matching tree entries are declared in BrowserView,
    # so the translation context requires this.
    # Merging this into BrowserView (and making the former more dynamic)
    # and using function pointers would go a long way towards
    # making the code more modular an easier to navigate thna this mess.

    return QCoreApplication.translate('BrowserView', text)

#-------------------------------------------------------------------------------
# End of Toolbox
#-------------------------------------------------------------------------------
