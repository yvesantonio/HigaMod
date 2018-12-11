# state file generated using paraview version 5.3.0

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# Create a new 'Render View'
renderView2 = CreateView('RenderView')
renderView2.ViewSize = [891, 590]
renderView2.AnnotationColor = [0.0, 0.0, 0.0]
renderView2.AxesGrid = 'GridAxes3DActor'
renderView2.OrientationAxesVisibility = 0
renderView2.OrientationAxesLabelColor = [0.0, 0.0, 0.0]
renderView2.OrientationAxesOutlineColor = [0.0, 0.0, 0.0]
renderView2.CenterOfRotation = [-2.76751066508041e-10, -2.76644485097677e-10, 4.12011900544167]
renderView2.StereoType = 0
renderView2.CameraPosition = [4.157579630668605, -18.211873706366443, 24.13439356533021]
renderView2.CameraFocalPoint = [-0.30739475390895354, 0.9622811494967001, 5.165149442790441]
renderView2.CameraViewUp = [-0.09033300072204026, 0.6897129972229421, 0.7184260090241011]
renderView2.CameraParallelScale = 8.56175748387977
renderView2.Background = [1.0, 1.0, 1.0]

# init the 'GridAxes3DActor' selected for 'AxesGrid'
renderView2.AxesGrid.XTitleColor = [0.0, 0.0, 0.0]
renderView2.AxesGrid.YTitleColor = [0.0, 0.0, 0.0]
renderView2.AxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
renderView2.AxesGrid.GridColor = [0.0, 0.0, 0.0]
renderView2.AxesGrid.XLabelColor = [0.0, 0.0, 0.0]
renderView2.AxesGrid.YLabelColor = [0.0, 0.0, 0.0]
renderView2.AxesGrid.ZLabelColor = [0.0, 0.0, 0.0]

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'Legacy VTK Reader'
fluidvtk = LegacyVTKReader(FileNames=['C:/Users/Simone/Desktop/AddBoundaryFlag/fluid.vtk'])

# create a new 'Extract Surface'
extractSurface1 = ExtractSurface(Input=fluidvtk)

# create a new 'SplineSource'
splineSource1 = SplineSource()
splineSource1.ParametricFunction = 'Spline'

# init the 'Spline' selected for 'ParametricFunction'
splineSource1.ParametricFunction.Points = [-2.34356558298757, 1.26036550406309, 9.23, -1.1445810552756, 0.327686343342759, 9.23, 0.564185815291419, -0.665820218279141, 9.23, 2.48538786452623, -1.57620544677045, 9.23, 3.20781246015624, -0.27037320869703, 9.23, 3.48477235923616, 1.57399160379305, 9.23, 2.41335500106081, 3.06389412892264, 9.23, 0.644993101043665, 3.67573482318455, 9.23, -0.563787905071035, 3.45439032889505, 9.23, -1.44055908632551, 2.8864907402776, 9.23]
splineSource1.ParametricFunction.Closed = 1

# create a new 'Delaunay 2D'
delaunay2D1 = Delaunay2D(Input=splineSource1)

# create a new 'Slice'
slice1 = Slice(Input=extractSurface1)
slice1.SliceType = 'Sphere'
slice1.SliceOffsetValues = [0.0]

# init the 'Sphere' selected for 'SliceType'
slice1.SliceType.Center = [-1.20696719887481, -2.24606757421074, 8.7708901850179]
slice1.SliceType.Radius = 1.8

# create a new 'Delaunay 2D'
delaunay2D2 = Delaunay2D(Input=slice1)

# ----------------------------------------------------------------
# setup the visualization in view 'renderView2'
# ----------------------------------------------------------------

# show data from splineSource1
splineSource1Display = Show(splineSource1, renderView2)
# trace defaults for the display properties.
splineSource1Display.Representation = 'Surface'
splineSource1Display.AmbientColor = [0.0, 0.0, 0.0]
splineSource1Display.ColorArrayName = [None, '']
splineSource1Display.OSPRayScaleFunction = 'PiecewiseFunction'
splineSource1Display.SelectOrientationVectors = 'None'
splineSource1Display.ScaleFactor = 0.588447046279907
splineSource1Display.SelectScaleArray = 'None'
splineSource1Display.GlyphType = 'Arrow'
splineSource1Display.PolarAxes = 'PolarAxesRepresentation'
splineSource1Display.GaussianRadius = 0.294223523139954
splineSource1Display.SetScaleArray = [None, '']
splineSource1Display.ScaleTransferFunction = 'PiecewiseFunction'
splineSource1Display.OpacityArray = [None, '']
splineSource1Display.OpacityTransferFunction = 'PiecewiseFunction'

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
splineSource1Display.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
splineSource1Display.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
splineSource1Display.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
splineSource1Display.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]

# show data from delaunay2D1
delaunay2D1Display = Show(delaunay2D1, renderView2)
# trace defaults for the display properties.
delaunay2D1Display.Representation = 'Surface With Edges'
delaunay2D1Display.AmbientColor = [0.0, 0.0, 0.0]
delaunay2D1Display.ColorArrayName = [None, '']
delaunay2D1Display.OSPRayScaleFunction = 'PiecewiseFunction'
delaunay2D1Display.SelectOrientationVectors = 'None'
delaunay2D1Display.ScaleFactor = 0.588447046279907
delaunay2D1Display.SelectScaleArray = 'None'
delaunay2D1Display.GlyphType = 'Arrow'
delaunay2D1Display.PolarAxes = 'PolarAxesRepresentation'
delaunay2D1Display.GaussianRadius = 0.294223523139954
delaunay2D1Display.SetScaleArray = [None, '']
delaunay2D1Display.ScaleTransferFunction = 'PiecewiseFunction'
delaunay2D1Display.OpacityArray = [None, '']
delaunay2D1Display.OpacityTransferFunction = 'PiecewiseFunction'

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
delaunay2D1Display.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
delaunay2D1Display.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
delaunay2D1Display.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
delaunay2D1Display.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]

# show data from delaunay2D2
delaunay2D2Display = Show(delaunay2D2, renderView2)
# trace defaults for the display properties.
delaunay2D2Display.Representation = 'Surface With Edges'
delaunay2D2Display.AmbientColor = [0.0, 0.0, 0.0]
delaunay2D2Display.ColorArrayName = [None, '']
delaunay2D2Display.OSPRayScaleFunction = 'PiecewiseFunction'
delaunay2D2Display.SelectOrientationVectors = 'None'
delaunay2D2Display.ScaleFactor = 0.346487258777642
delaunay2D2Display.SelectScaleArray = 'None'
delaunay2D2Display.GlyphType = 'Arrow'
delaunay2D2Display.PolarAxes = 'PolarAxesRepresentation'
delaunay2D2Display.GaussianRadius = 0.173243629388821
delaunay2D2Display.SetScaleArray = [None, '']
delaunay2D2Display.ScaleTransferFunction = 'PiecewiseFunction'
delaunay2D2Display.OpacityArray = [None, '']
delaunay2D2Display.OpacityTransferFunction = 'PiecewiseFunction'

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
delaunay2D2Display.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
delaunay2D2Display.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
delaunay2D2Display.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
delaunay2D2Display.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]

# show data from fluidvtk
fluidvtkDisplay = Show(fluidvtk, renderView2)
# trace defaults for the display properties.
fluidvtkDisplay.Representation = 'Surface With Edges'
fluidvtkDisplay.AmbientColor = [0.0, 0.0, 0.0]
fluidvtkDisplay.ColorArrayName = [None, '']
fluidvtkDisplay.Opacity = 0.23
fluidvtkDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
fluidvtkDisplay.SelectOrientationVectors = 'None'
fluidvtkDisplay.ScaleFactor = 1.0219761073589326
fluidvtkDisplay.SelectScaleArray = 'None'
fluidvtkDisplay.GlyphType = 'Arrow'
fluidvtkDisplay.PolarAxes = 'PolarAxesRepresentation'
fluidvtkDisplay.ScalarOpacityUnitDistance = 0.5270777028912511
fluidvtkDisplay.GaussianRadius = 0.5109880536794663
fluidvtkDisplay.SetScaleArray = [None, '']
fluidvtkDisplay.ScaleTransferFunction = 'PiecewiseFunction'
fluidvtkDisplay.OpacityArray = [None, '']
fluidvtkDisplay.OpacityTransferFunction = 'PiecewiseFunction'

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
fluidvtkDisplay.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
fluidvtkDisplay.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
fluidvtkDisplay.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
fluidvtkDisplay.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]

# show data from extractSurface1
extractSurface1Display = Show(extractSurface1, renderView2)
# trace defaults for the display properties.
extractSurface1Display.Representation = 'Surface'
extractSurface1Display.AmbientColor = [0.0, 0.0, 0.0]
extractSurface1Display.ColorArrayName = [None, '']
extractSurface1Display.OSPRayScaleFunction = 'PiecewiseFunction'
extractSurface1Display.SelectOrientationVectors = 'None'
extractSurface1Display.ScaleFactor = 1.0219761073589326
extractSurface1Display.SelectScaleArray = 'None'
extractSurface1Display.GlyphType = 'Arrow'
extractSurface1Display.PolarAxes = 'PolarAxesRepresentation'
extractSurface1Display.GaussianRadius = 0.5109880536794663
extractSurface1Display.SetScaleArray = [None, '']
extractSurface1Display.ScaleTransferFunction = 'PiecewiseFunction'
extractSurface1Display.OpacityArray = [None, '']
extractSurface1Display.OpacityTransferFunction = 'PiecewiseFunction'

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
extractSurface1Display.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
extractSurface1Display.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
extractSurface1Display.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
extractSurface1Display.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]

# show data from slice1
slice1Display = Show(slice1, renderView2)
# trace defaults for the display properties.
slice1Display.Representation = 'Surface'
slice1Display.AmbientColor = [0.0, 0.0, 0.0]
slice1Display.ColorArrayName = [None, '']
slice1Display.OSPRayScaleFunction = 'PiecewiseFunction'
slice1Display.SelectOrientationVectors = 'None'
slice1Display.ScaleFactor = 0.3464872587776413
slice1Display.SelectScaleArray = 'None'
slice1Display.GlyphType = 'Arrow'
slice1Display.PolarAxes = 'PolarAxesRepresentation'
slice1Display.GaussianRadius = 0.17324362938882065
slice1Display.SetScaleArray = [None, '']
slice1Display.ScaleTransferFunction = 'PiecewiseFunction'
slice1Display.OpacityArray = [None, '']
slice1Display.OpacityTransferFunction = 'PiecewiseFunction'

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
slice1Display.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
slice1Display.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
slice1Display.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
slice1Display.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]

# ----------------------------------------------------------------
# finally, restore active source
SetActiveSource(fluidvtk)
# ----------------------------------------------------------------
