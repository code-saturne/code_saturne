#------------------------------------------------------------------------------
#   This file is part of the Code_Saturne Kernel, element of the
#   Code_Saturne CFD tool.
#
#   Copyright (C) 2009 EDF S.A., France
#
#   The Code_Saturne Kernel is free software; you can redistribute it
#   and/or modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2 of
#   the License, or (at your option) any later version.
#
#   The Code_Saturne Kernel is distributed in the hope that it will be
#   useful, but WITHOUT ANY WARRANTY; without even the implied warranty
#   of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public Licence
#   along with the Code_Saturne Kernel; if not, write to the
#   Free Software Foundation, Inc.,
#   51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
#-------------------------------------------------------------------------------

all-local:
	@$(PYTHON) $(top_srcdir)/gui/sbin/mkpyqt.py --build $(srcdir)

clean-local:
	@$(PYTHON) $(top_srcdir)/gui/sbin/mkpyqt.py --clean $(srcdir)

# Install Python library

pagesdir = $(pkgpythondir)/Pages
pages_PYTHON = \
__init__.py \
resources_pages_rc.py \
AnalysisFeaturesForm.py \
AnalysisFeaturesView.py \
AtmosphericFlowsForm.py \
AtmosphericFlowsModel.py \
AtmosphericFlowsView.py \
BatchRunningAdvancedOptionsDialogForm.py \
BatchRunningForm.py \
BatchRunningModel.py \
BatchRunningPBSJobManagementDialogForm.py \
BatchRunningUserFilesDialogForm.py \
BatchRunningView.py \
BodyForcesForm.py \
BodyForcesModel.py \
BodyForcesView.py \
BoundaryConditionsCoalInletForm.py \
BoundaryConditionsCoalInletView.py \
BoundaryConditionsForm.py \
BoundaryConditionsMeteoForm.py \
BoundaryConditionsMeteoView.py \
BoundaryConditionsMobileMeshForm.py \
BoundaryConditionsMobileMeshView.py \
BoundaryConditionsRoughWallForm.py \
BoundaryConditionsRoughWallView.py \
BoundaryConditionsScalarsForm.py \
BoundaryConditionsScalarsView.py \
BoundaryConditionsSlidingWallForm.py \
BoundaryConditionsSlidingWallView.py \
BoundaryConditionsTurbulenceInletForm.py \
BoundaryConditionsTurbulenceInletView.py \
BoundaryConditionsVelocityInletForm.py \
BoundaryConditionsVelocityInletView.py \
BoundaryConditionsView.py \
Boundary.py \
CoalCombustionForm.py \
CoalCombustionModel.py \
CoalCombustionView.py \
CoalThermoChemistry.py \
CommonCombustion.py \
CurrentSpeciesForm.py \
CurrentSpeciesModel.py \
CurrentSpeciesView.py \
DefineUserScalarsForm.py \
DefineUserScalarsModel.py \
DefineUserScalarsView.py \
ElectricalModelsModel.py \
FacesSelectionForm.py \
FacesSelectionView.py \
FluidCharacteristicsForm.py \
FluidCharacteristicsModel.py \
FluidCharacteristicsView.py \
FluidStructureInteractionAdvancedOptionsDialogForm.py \
FluidStructureInteractionForm.py \
FluidStructureInteractionModel.py \
FluidStructureInteractionView.py \
GasCombustionModel.py \
HeadLossesForm.py \
HeadLossesModel.py \
HeadLossesView.py \
IdentityAndPathesForm.py \
IdentityAndPathesModel.py \
IdentityAndPathesView.py \
InitializationForm.py \
InitializationModel.py \
InitializationView.py \
LagrangianAdvancedOptionsDialogForm.py \
LagrangianBoundariesForm.py \
LagrangianBoundariesModel.py \
LagrangianBoundariesView.py \
LagrangianForm.py \
LagrangianModel.py \
LagrangianOutputForm.py \
LagrangianOutputModel.py \
LagrangianOutputView.py \
LagrangianStatisticsForm.py \
LagrangianStatisticsModel.py \
LagrangianStatisticsView.py \
LagrangianView.py \
LocalizationForm.py \
LocalizationModel.py \
LocalizationView.py \
MatisseCustomForm.py \
MatisseCustomView.py \
MatisseGeomForm.py \
MatisseGeomModel.py \
MatisseGeomView.py \
MatisseHydrauForm.py \
MatisseHydrauModel.py \
MatisseHydrauView.py \
MatisseModel.py \
MatisseNetworkForm.py \
MatisseNetworkModel.py \
MatisseNetworkView.py \
MatisseRangeDescriptionForm.py \
MatisseRangeDescriptionModel.py \
MatisseRangeDescriptionView.py \
MatisseThermicForm.py \
MatisseThermicModel.py \
MatisseThermicView.py \
MatisseTypeForm.py \
MatisseTypeModel.py \
MatisseTypeView.py \
MemoryAllocationForm.py \
MemoryAllocationModel.py \
MemoryAllocationView.py \
MeshQualityCriteriaLogDialogForm.py \
MobileMeshForm.py \
MobileMeshModel.py \
MobileMeshView.py \
NumericalParamEquationForm.py \
NumericalParamEquationModel.py \
NumericalParamEquationView.py \
NumericalParamGlobalForm.py \
NumericalParamGlobalModel.py \
NumericalParamGlobalView.py \
OutputControlForm.py \
OutputControlModel.py \
OutputControlView.py \
OutputSurfacicVariablesForm.py \
OutputSurfacicVariablesModel.py \
OutputSurfacicVariablesView.py \
OutputVolumicVariablesForm.py \
OutputVolumicVariablesModel.py \
OutputVolumicVariablesView.py \
PreProcessingInformationsView.py \
ProfilesForm.py \
ProfilesModel.py \
ProfilesView.py \
QMeiEditorForm.py \
QMeiEditorView.py \
RadiativeBoundariesForm.py \
RadiativeBoundariesView.py \
ReferenceValuesForm.py \
ReferenceValuesModel.py \
ReferenceValuesView.py \
SalomeHandler.py \
SolutionDomainForm.py \
SolutionDomainMeshFormatDialogForm.py \
SolutionDomainModel.py \
SolutionDomainView.py \
SolutionVerifForm.py \
SolutionVerifView.py \
StartRestartAdvancedDialogForm.py \
StartRestartForm.py \
StartRestartModel.py \
StartRestartView.py \
SteadyManagementForm.py \
SteadyManagementModel.py \
SteadyManagementView.py \
SyrthesForm.py \
SyrthesView.py \
ThermalRadiationAdvancedDialogForm.py \
ThermalRadiationForm.py \
ThermalRadiationModel.py \
ThermalRadiationView.py \
ThermalScalarForm.py \
ThermalScalarModel.py \
ThermalScalarView.py \
TimeAveragesForm.py \
TimeAveragesModel.py \
TimeAveragesView.py \
TimeStepForm.py \
TimeStepModel.py \
TimeStepView.py \
TurbulenceAdvancedOptionsDialogForm.py \
TurbulenceForm.py \
TurbulenceModel.py \
TurbulenceView.py \
UserArraysForm.py \
UserArraysModel.py \
UserArraysView.py \
UserScalarPropertiesForm.py \
UserScalarPropertiesView.py \
VerifyExistenceLabelDialogForm.py \
VerifyExistenceLabelDialogView.py \
WelcomeForm.py \
WelcomeView.py
