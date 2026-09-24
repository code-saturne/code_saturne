# script-version: 2.0
# Catalyst state generated using paraview version 6.1.1-916-g0228d803a6

# Generated and usable with code_saturne 12_Tee_junction/CASE3 tutorial

import paraview
paraview.compatibility.major = 6
paraview.compatibility.minor = 1

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.Set(
    ViewSize=[1021, 548],
    CenterOfRotation=[-3.3527612686157227e-08, 0.4749999940395355, 0.08750000223517418],
    CameraPosition=[0.6309863357200569, -0.8384419310732949, 0.24934816622958506],
    CameraFocalPoint=[-0.1500084685659077, 0.3495984185819664, -0.027471373185893894],
    CameraViewUp=[-0.20926453945165482, 0.08932709104977768, 0.9737705188249806],
)

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'EnSight Reader'
catalyst = EnSightReader(registrationName='catalyst', CaseFileName='postprocessing/RESULTS_FLUID_DOMAIN.case')
catalyst.CellArrays = ['Velocity', 'Pressure', 'k', 'epsilon', 'TempC', 'TurbVisc', 'CourantNb', 'FourierNb', 'total_pressure', 'Local_Time_Step', 'mpi_rank_id']

# create a new 'Ghost Cells'
ghostCells1 = GhostCells(registrationName='GhostCells1', Input=catalyst)

# create a new 'Cell Data to Point Data'
cellDatatoPointData1 = CellDatatoPointData(registrationName='CellDatatoPointData1', Input=ghostCells1)

# create a new 'Stream Tracer'
streamTracer1 = StreamTracer(registrationName='StreamTracer1', Input=cellDatatoPointData1,
    SeedType='Point Cloud')
streamTracer1.Set(
    Vectors=['POINTS', 'Velocity'],
    MaximumStreamlineLength=1.550000011920929,
    ComputeVorticity=0,
)

# init the 'Point Cloud' selected for 'SeedType'
streamTracer1.SeedType.Set(
    Center=[0.0, -0.3, 0.0],
    Radius=0.1550000011920929,
)

# create a new 'Stream Tracer'
streamTracer2 = StreamTracer(registrationName='StreamTracer2', Input=cellDatatoPointData1,
    SeedType='Point Cloud')
streamTracer2.Set(
    Vectors=['POINTS', 'Velocity'],
    MaximumStreamlineLength=1.550000011920929,
)

# init the 'Point Cloud' selected for 'SeedType'
streamTracer2.SeedType.Set(
    Center=[0.05, 0.0, 0.25],
    NumberOfPoints=20,
    Radius=0.05,
)

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from streamTracer1
streamTracer1Display = Show(streamTracer1, renderView1, 'GeometryRepresentation')

# get color transfer function/color map for 'Velocity'
velocityLUT = GetColorTransferFunction('Velocity')
velocityLUT.Set(
    RGBPoints=GenerateRGBPoints(
        range_min=0.3324082844862145,
        range_max=2.140373181426344,
    ),
    ScalarRangeInitialized=1.0,
)

# trace defaults for the display properties.
streamTracer1Display.Set(
    Representation='Surface',
    ColorArrayName=['POINTS', 'Velocity'],
    LookupTable=velocityLUT,
)

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
streamTracer1Display.ScaleTransferFunction.Points = [-198.663818359375, 0.0, 0.5, 0.0, 273.2964782714844, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
streamTracer1Display.OpacityTransferFunction.Points = [-198.663818359375, 0.0, 0.5, 0.0, 273.2964782714844, 1.0, 0.5, 0.0]

# show data from streamTracer2
streamTracer2Display = Show(streamTracer2, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
streamTracer2Display.Set(
    Representation='Surface',
    ColorArrayName=['POINTS', 'Velocity'],
    LookupTable=velocityLUT,
)

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
streamTracer2Display.ScaleTransferFunction.Points = [-94.8069839477539, 0.0, 0.5, 0.0, 188.6370391845703, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
streamTracer2Display.OpacityTransferFunction.Points = [-94.8069839477539, 0.0, 0.5, 0.0, 188.6370391845703, 1.0, 0.5, 0.0]

# show data from ghostCells1
ghostCells1Display = Show(ghostCells1, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
ghostCells1Display.Set(
    Representation='Surface',
    ColorArrayName=[None, ''],
    Opacity=0.1,
    Assembly='Hierarchy',
)

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
ghostCells1Display.ScaleTransferFunction.Points = [-1399.7635498046875, 0.0, 0.5, 0.0, 317.8004150390625, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
ghostCells1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]
ghostCells1Display.OpacityTransferFunction.Points = [-1399.7635498046875, 0.0, 0.5, 0.0, 317.8004150390625, 1.0, 0.5, 0.0]

# setup the color legend parameters for each legend in this view

# get color legend/bar for velocityLUT in view renderView1
velocityLUTColorBar = GetScalarBar(velocityLUT, renderView1)
velocityLUTColorBar.Set(
    AutoOrient=0,
    Orientation='Horizontal',
    WindowLocation='Any Location',
    Position=[0.55162521075924, 0.07481751824817519],
    Title='Velocity',
    ComponentTitle='Magnitude',
    ScalarBarLength=0.1858394160583945,
)

# set color bar visibility
velocityLUTColorBar.Visibility = 1

# show color legend
streamTracer1Display.SetScalarBarVisibility(renderView1, True)

# show color legend
streamTracer2Display.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# setup color maps and opacity maps used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get opacity transfer function/opacity map for 'Velocity'
velocityPWF = GetOpacityTransferFunction('Velocity')
velocityPWF.Set(
    Points=[0.3324082844862145, 0.0, 0.5, 0.0, 2.140373181426344, 1.0, 0.5, 0.0],
    ScalarRangeInitialized=1,
)

# ----------------------------------------------------------------
# setup extractors
# ----------------------------------------------------------------

# create extractor
pNG1 = CreateExtractor('PNG', renderView1, registrationName='PNG1')
# trace defaults for the extractor.
# init the 'PNG' selected for 'Writer'
pNG1.Writer.Set(
    FileName='stream_lines1_{timestep:06d}.png',
    ImageResolution=[1021, 548],
    OverrideColorPalette='WhiteBackground',
    Format='PNG',
)

# ----------------------------------------------------------------
# setup camera links
# ----------------------------------------------------------------

# ------------------------------------------------------------------------------
# Catalyst options
from paraview import catalyst
options = catalyst.Options()

# ------------------------------------------------------------------------------
if __name__ == '__main__':
    from paraview.simple import SaveExtractsUsingCatalystOptions
    # Code for non in-situ environments; if executing in post-processing
    # i.e. non-Catalyst mode, let's generate extracts using Catalyst options
    SaveExtractsUsingCatalystOptions(options)
