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
renderView2.ViewSize = [815, 590]
renderView2.AnnotationColor = [0.0, 0.0, 0.0]
renderView2.AxesGrid = 'GridAxes3DActor'
renderView2.OrientationAxesLabelColor = [0.0, 0.0, 0.0]
renderView2.OrientationAxesOutlineColor = [0.0, 0.0, 0.0]
renderView2.CenterOfRotation = [0.5807347297668457, 1.0280951261520386, 7.729999542236328]
renderView2.StereoType = 0
renderView2.CameraPosition = [7.130284947474139, -8.044318184568423, 24.05293866563277]
renderView2.CameraFocalPoint = [0.5807347297668279, 1.0280951261520503, 7.72999954223633]
renderView2.CameraViewUp = [-0.4956998922762894, 0.6592963177874992, 0.5653405895114201]
renderView2.CameraParallelScale = 4.233079213673913
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
legacyVTKReader1 = LegacyVTKReader(FileNames=['C:/Users/Simone/Desktop/AddBoundaryFlag/fluid.vtk'])

# create a new 'Extract Surface'
extractSurface1 = ExtractSurface(Input=legacyVTKReader1)

# create a new 'SplineSource'
splineSource1 = SplineSource()
splineSource1.ParametricFunction = 'Spline'

# init the 'Spline' selected for 'ParametricFunction'
splineSource1.ParametricFunction.Points = [-2.34356558298757, 1.26036550406309, 9.23, -1.1445810552756, 0.327686343342759, 9.23, 0.564185815291419, -0.665820218279141, 9.23, 2.48538786452623, -1.57620544677045, 9.23, 3.20781246015624, -0.27037320869703, 9.23, 3.48477235923616, 1.57399160379305, 9.23, 2.41335500106081, 3.06389412892264, 9.23, 0.644993101043665, 3.67573482318455, 9.23, -0.563787905071035, 3.45439032889505, 9.23, -1.44055908632551, 2.8864907402776, 9.23]
splineSource1.ParametricFunction.Closed = 1

# create a new 'Linear Extrusion'
linearExtrusion1 = LinearExtrusion(Input=splineSource1)
linearExtrusion1.ScaleFactor = 3.0
linearExtrusion1.Vector = [0.0, 0.0, -1.0]

# create a new 'Extract Surface'
extractSurface2 = ExtractSurface(Input=linearExtrusion1)

# create a new 'Triangulate'
triangulate1 = Triangulate(Input=extractSurface2)

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

# show data from triangulate1
triangulate1Display = Show(triangulate1, renderView2)
# trace defaults for the display properties.
triangulate1Display.Representation = 'Wireframe'
triangulate1Display.AmbientColor = [0.0, 0.0, 0.0]
triangulate1Display.ColorArrayName = [None, '']
triangulate1Display.OSPRayScaleFunction = 'PiecewiseFunction'
triangulate1Display.SelectOrientationVectors = 'None'
triangulate1Display.ScaleFactor = 0.5884470462799073
triangulate1Display.SelectScaleArray = 'None'
triangulate1Display.GlyphType = 'Arrow'
triangulate1Display.PolarAxes = 'PolarAxesRepresentation'
triangulate1Display.GaussianRadius = 0.29422352313995365
triangulate1Display.SetScaleArray = [None, '']
triangulate1Display.ScaleTransferFunction = 'PiecewiseFunction'
triangulate1Display.OpacityArray = [None, '']
triangulate1Display.OpacityTransferFunction = 'PiecewiseFunction'

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
triangulate1Display.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
triangulate1Display.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
triangulate1Display.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
triangulate1Display.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]

# ----------------------------------------------------------------
# finally, restore active source
SetActiveSource(triangulate1)
# ----------------------------------------------------------------
