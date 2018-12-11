#!/usr/bin/env python

## Date:      $Date: 2014/02/18 12:35:13 $

## Note: this class was contributed by 
##       Elena Faggiano (elena.faggiano@gmail.com)
##       Universita di Pavia

##IMPORTANT: If you are going to publish with the help of this script please contact Elena Faggiano for a publication to cite.

#vmtk installation required (http://www.vmtk.org/documentation/installation.html)

from __future__ import absolute_import #NEEDS TO STAY AS TOP LEVEL MODULE FOR Py2-3 COMPATIBILITY
import math
import sys
import vtk

from vmtk import pypes
from vmtk import vmtkscripts
from vmtk import vtkvmtk
#from vmtk import vmtkrenderer


fsimeshgenerator = 'FSImeshgenerator'

class FSImeshgenerator(pypes.pypeScript):

    def __init__(self):

        pypes.pypeScript.__init__(self)

        self.Surface = None
        self.Centerlines = None
        self.SkipCapAndRemesh = 0
        self.RemeshCapsOnly = 0
        self.InterfaceLabel = 1
        self.LifeVFluidMeshName = None
        self.LifeVSolidMeshName = None
        self.VTKPrefixFilesName =None
        self.ScaleFactor = 1

        #TOGGLE INTERACTIVE
        self.SurfaceInteractive = 1
        self.FluidInteractive = 1
        self.SolidInteractive = 1

        #SURFACE
        self.RemeshingMethod = 'radius' #'constant' , 'edgelengtharray'
        self.TargetEdgeLengthArrayName = 'TargetEdgelength'
        self.TargetEdgeLength = 1.0
        self.TargetAlphaFactor = 0.2
        self.TargetBetaFactor = 1.0

        #FLUID
        self.Fluid = 1
        self.VolumeElementScaleFactor = 0.8

        #BOUNDARY LAYER
        self.BoundaryLayer = 0
        self.BoundaryLayerNumberOfSubLayers = 2
        self.SubLayerRatio = 0.8
        self.BoundaryLayerNumberOfSubsteps = 2000
        self.BoundaryLayerMethod = 'radius' #o radius
        self.BoundaryLayerConstantThickness = 0.1
        self.BoundaryLayerAlphaFactor = 0.1
        self.BoundaryLayerBetaFactor = 0.5

        #CAPS (NEEDED ONLY IF SkipCapAndRemesh == 1 AND BoundaryLayer == 0)
        self.CapRemeshingMethod = 'radius'
        self.CapTargetEdgeLength = 1.0
        self.CapTargetAlphaFactor = 0.2
        self.CapTargetBetaFactor = 1.0

        #SOLID
        self.Solid = 1
        self.SolidNumberOfSubLayers = 2
        self.SolidMethod = 'radius' #o radius
        self.SolidNumberOfSubsteps = 2000
        self.SolidConstantThickness = 1.0
        self.SolidAlphaFactor = 0.2
        self.NormalsCorrected = 0

        #OUTPUTS
        self.RemeshedSurface = None
        self.FluidMesh = None
        self.SolidMesh = None
        self.ComputedCenterline = None


        self.SetScriptName('LifeVFSIgeneratorOpt')
        self.SetScriptDoc('the code generates a fluid mesh and a structure mesh in LifeV format starting from an interface mesh. Multiple areas delineation is also supported.')
        self.SetInputMembers([
            ['Surface','i','vtkPolyData',1,'','the input surface','vmtksurfacereader'],
            ['Centerlines', 'centerlines','vtkPolyData',1,'','centerlines previously calculated (optional)','vmtksurfacereader'],
            ['SkipCapAndRemesh','skipcapremesh','bool',1,'','set 1 if your input is a capped and remeshed surface with labels'],
            ['RemeshCapsOnly','remeshcapsonly','bool',1,'','set 1 if you want to remesh only caps'],
            ['InterfaceLabel','ilabel','float',1,'(0.0,)','label of the interface (needed if SkipCapAndRemesh = 1)'],
            ['LifeVFluidMeshName','fluidname','str',1,'','fluid mesh file name (with .mesh extension)'],
            ['LifeVSolidMeshName','solidname','str',1,'','solid mesh file name (with .mesh extension)'],
            ['VTKPrefixFilesName','vtkprefixname','str',1,'','vtk files prefix : OutputName (without extention) for the generation of the files: OutputName_centerlines.vtp OutputName_interface.vtp OutputName_labeledinterface.vtp OutputName_fluid.vtu and OutputName_solid.vtu'],
            ['ScaleFactor','scale','float',1,'(0.0,)','scaling parameter to scale the mesh coordinates (e.g. from mm to cm -scale 0.1)'],           
            ['SurfaceInteractive','surfaceinteractive','bool',1,'','toggle interactive generation of the surface'],
            ['FluidInteractive','fluidinteractive','bool',1,'','toggle interactive generation of the fluid'],
            ['SolidInteractive','solidinteractive','bool',1,'','toggle interactive generation of the solid'],
            ['Fluid','fluid','bool',1,'','toggle fluid generation'],
            ['RemeshingMethod','remeshingmethod','str',1,'["constant","radius", "edgelengtharray"]','choose between constant or radius dependent edgelength (surface remeshing step)'],
            ['TargetEdgeLengthArrayName','edgelengtharray','str',1],
            ['TargetEdgeLength','edgelength','float',1,'(0.0,)','target constant edgelength (surface remeshing step)'],
            ['TargetAlphaFactor','alpha','float',1,'(0.0,)','target alpha factor (surface remeshing step with h = aplha * radius ^ beta)'],
            ['TargetBetaFactor','beta','float',1,'(0.0,)','target beta factor (surface remeshing step with h = aplha * radius ^ beta)'],
            ['VolumeElementScaleFactor','volumeelementfactor','float',1,'(0.0,)','target volume element scale factor'],
            ['BoundaryLayer','boundarylayer','bool',1,'','toggle boundary layer generation'],
            ['BoundaryLayerNumberOfSubLayers','blsublayers','int',1,'(0,)','number of sublayers in boundary layer'],
            ['BoundaryLayerNumberOfSubsteps','blsubsteps','int',1,'(0,)','number of substeps for smoothly propagating the boundary layer'],
            ['SubLayerRatio','blsublayerratio','float',1,'(0.0,)','ratio between succesive layers'],
            ['BoundaryLayerMethod','blmethod','str',1,'["constant","radius"]','choose between constant or radius dependent boundary layer thickess'],
            ['BoundaryLayerConstantThickness','blconstantthickness','float',1,'(0.0,)','target constant thickness'],
            ['BoundaryLayerAlphaFactor','blalphafactor','float',1,'(0.0,)','target alpha factor (thickess t = aplha * radius ^ beta)'],
            ['BoundaryLayerBetaFactor','blbetafactor','float',1,'(0.0,)','target beta factor (thickess t = aplha * radius ^ beta)'],
            ['CapRemeshingMethod','capremeshingmethod','str',1,'["constant","radius"]','choose between constant or radius dependent edgelength for the cap generation (if skipremeshing = 1 and boundarylayer = 1)'],
            ['CapTargetEdgeLength','capedgelength','float',1,'(0.0,)','target constant edgelength for the caps (if skipcapremesh = 1 and boundarylayer = 1)'],
            ['CapTargetAlphaFactor','capalpha','float',1,'(0.0,)','target alpha factor (if skipcapremesh = 1 and boundarylayer = 1)'],
            ['CapTargetBetaFactor','capbeta','float',1,'(0.0,)','target beta factor (if skipcapremesh = 1 and boundarylayer = 1)'],
            ['Solid','solid','bool',1,'','toggle solid mesh generation'],
            ['SolidNumberOfSubLayers','solidsublayers','int',1,'(0,)','number of layer in solid'],
            ['SolidNumberOfSubsteps','solidsubsteps','int',1,'(0,)','number of substeps for smoothly propagating the boundary layer'],
            ['SolidMethod','solidmethod','str',1,'["constant","radius"]','choose between constant or radius dependent solid thickess'],
            ['SolidConstantThickness','solidconstantthickness','float',1,'(0.0,)','target constant thickess'],
            ['SolidAlphaFactor','solidalphafactor','float',1,'(0.0,)','target alpha factor (thickess t = aplha * radius )'],
            ['NormalsCorrected','normalscorrected','bool',1,'','toggle use of normals on rings perpendicular to inlet/outlet in the generation of solid'],
            ])
        self.SetOutputMembers([
            ['RemeshedSurface','remeshedsurface','vtkPolyData',1,'','the output surface','vmtksurfacewriter'],
            ['ComputedCenterline','ccenterline','vtkPolyData',1,'','the computed centerline','vmtksurfacewriter'],
            ['FluidMesh','fluidmesh','vtkUnstructuredGrid',1,'','the output fluid mesh','vmtkmeshwriter'],
            ['SolidMesh','solidmesh','vtkUnstructuredGrid',1,'','the output solid mesh','vmtkmeshwriter'],
            ])

    def WritePolyData(self,surface,name):
        writer = vtk.vtkXMLPolyDataWriter()
        print 'Writing surface:', name, '\n'
        writer.SetFileName(name)
        writer.SetInputData(surface)
        writer.Update()

    def WriteUnstructuredGrid(self,mesh,name):
        writer = vtk.vtkXMLUnstructuredGridWriter()
        print 'Writing mesh:', name, '\n'
        writer.SetFileName(name)
        writer.SetInputData(mesh)
        writer.Update()

    def ComputeSmartRadius(self,surfaceIn,centerlines,cellentityidsarrayname,radiusarrayname,distancetocenterlinearrayname,smartdistancearrayname):
        distanceToCenterlinesFilter = vtkvmtk.vtkvmtkPolyDataDistanceToCenterlines()
        distanceToCenterlinesFilter.SetInputData(surfaceIn)
        distanceToCenterlinesFilter.SetCenterlines(centerlines)
        distanceToCenterlinesFilter.SetUseRadiusInformation(1)            
        distanceToCenterlinesFilter.SetEvaluateCenterlineRadius(1)
        distanceToCenterlinesFilter.SetCenterlineRadiusArrayName(radiusarrayname)
        distanceToCenterlinesFilter.SetDistanceToCenterlinesArrayName(distancetocenterlinearrayname)
        distanceToCenterlinesFilter.Update()
        surface = distanceToCenterlinesFilter.GetOutput()
        distanceArray = vtk.vtkFloatArray()
        distanceArray.SetName(smartdistancearrayname)
        distanceArray.SetNumberOfComponents(1)
        distanceArray.SetNumberOfTuples(surface.GetNumberOfPoints())
        surface.GetPointData().AddArray(distanceArray)                          
        centerlineArray = surface.GetPointData().GetArray(distancetocenterlinearrayname)
        radiusArray = surface.GetPointData().GetArray(radiusarrayname)
        for i in range (surface.GetNumberOfPoints()):
            centerlineval = centerlineArray.GetComponent(i,0)
            radius = radiusArray.GetComponent(i,0)
            if centerlineval > 1.4 * radius:
                distanceArray.SetTuple1(i,1.4*radius)
            elif centerlineval < 0.9 * radius:
                distanceArray.SetTuple1(i,radius)
            else:
                distanceArray.SetTuple1(i,centerlineval)
        return surface

    def SmoothPointDataSimple(self,surface,arrayname,connexity):
        array=surface.GetPointData().GetArray(arrayname)        
        extractEdges = vtk.vtkExtractEdges()
        extractEdges.SetInputData(surface)
        extractEdges.Update()        
        surfEdges = extractEdges.GetOutput()
        
        if connexity == 1:
            for i in range (surfEdges.GetNumberOfPoints()):
                cells = vtk.vtkIdList()
                surfEdges.GetPointCells(i,cells)
                vval = 0
                ddd = 0
                d = 0
                N = 0
                for j in range (cells.GetNumberOfIds()):
                    points = vtk.vtkIdList()
                    surfEdges.GetCellPoints(cells.GetId(j),points)
                    for k in range (points.GetNumberOfIds()):
                        if points.GetId(k) != i:
                            d = math.sqrt(vtk.vtkMath.Distance2BetweenPoints(surface.GetPoint(i),surface.GetPoint(points.GetId(k))))
                            dd = 1/d
                            val = array.GetComponent(points.GetId(k),0)#*dd
                            N = N+1
                            vval = vval + val
                            ddd = ddd + dd
                val = array.GetComponent(i,0)
                vval = vval + val           
                newval = vval / (N+1)
                array.SetTuple1(i,newval)            
        elif connexity == 2:
            for i in range (surfEdges.GetNumberOfPoints()):
                cells = vtk.vtkIdList()
                surfEdges.GetPointCells(i,cells)
                pointlist = vtk.vtkIdList()
                vval = 0
                ddd = 0
                d = 0
                N = 0
                for j in range (cells.GetNumberOfIds()):
                    points = vtk.vtkIdList()
                    surfEdges.GetCellPoints(cells.GetId(j),points)
                    for k in range (points.GetNumberOfIds()):
                        if points.GetId(k) != i:
                            pointlist.InsertUniqueId(points.GetId(k))
                            cells2 = vtk.vtkIdList()
                            surfEdges.GetPointCells(i,cells2)
                            for p in range (cells2.GetNumberOfIds()):
                                points2 = vtk.vtkIdList()
                                surfEdges.GetCellPoints(cells2.GetId(j),points2)
                                for q in range (points2.GetNumberOfIds()):
                                    if points2.GetId(k) != i:
                                        pointlist.InsertUniqueId(points2.GetId(k))                                              
                N = pointlist.GetNumberOfIds()
                for j in range (pointlist.GetNumberOfIds()):
                    d = math.sqrt(vtk.vtkMath.Distance2BetweenPoints(surface.GetPoint(i),surface.GetPoint(pointlist.GetId(j))))
                    dd = 1/d
                    val = array.GetComponent(pointlist.GetId(j),0)#*dd
                    vval = vval + val
                    ddd = ddd + dd
                val = array.GetComponent(i,0)
                vval = vval + val           
                newval = vval / (N+1)
                array.SetTuple1(i,newval)
        else:
            print "error, wrong connexity"
        return surface

    def Execute(self): 
        if (self.Surface == None):
            self.PrintError('Error: no Surface.') 

        if self.Fluid == 1:
            if (self.LifeVFluidMeshName == None):
                self.PrintError('Error: no Fluid mesh name.') 
        if self.Solid == 1:
            if (self.LifeVSolidMeshName == None):
                self.PrintError('Error: no Solid mesh name.')

        #--------- classic vmtk names ----------------
        CellEntityIdsArrayName = 'CellEntityIds'
        DistanceToCenterlinesArrayName = 'DistanceToCenterlines'
        RadiusArrayName = 'MaximumInscribedSphereRadius'
        SmartDistanceArrayName = 'SmartDistance'
        SizingFunctionArrayName = 'VolumeSizingFunction'
        #-------IDS PARAMETERS ----------
        #Mesh Ids
        InterfaceId = 200
        ExternalId = 210
        #--------------------------------
        triangleCellType = 5
        tetraCellType = 10
        #--------------------------------
        surface = self.Surface

        if (self.SkipCapAndRemesh == 0):
            meshing="n"
            if (self.Centerlines == None):
                centerlinesFilter = vmtkscripts.vmtkCenterlines()
                centerlinesFilter.Surface=surface
                centerlinesFilter.SeedSelectorName = 'openprofiles'
                centerlinesFilter.AppendEndPoints = 1
                centerlinesFilter.Execute()
                centerlines=centerlinesFilter.Centerlines
            else:
                centerlines = self.Centerlines
        else:
            ILabel = self.InterfaceLabel
            print "interface label", ILabel
            meshing = "y"
            print "interface mesh extraction\n"
            threshold = vtk.vtkThreshold()
            threshold.SetInputData(surface)
            threshold.ThresholdBetween(self.InterfaceLabel-0.5,self.InterfaceLabel+0.5)
            threshold.SetInputArrayToProcess(0, 0, 0, 1, CellEntityIdsArrayName)
            threshold.Update()
            meshToSurface = vmtkscripts.vmtkMeshToSurface()
            meshToSurface.Mesh = threshold.GetOutput()
            meshToSurface.Execute()
            interfaceSurface = meshToSurface.Surface
            surfaceRemeshed = surface
            if (self.Centerlines == None):
                centerlinesFilter = vmtkscripts.vmtkCenterlines()
                centerlinesFilter.Surface = interfaceSurface
                centerlinesFilter.SeedSelectorName = 'openprofiles'
                centerlinesFilter.AppendEndPoints = 1
                centerlinesFilter.Execute()
                centerlines = centerlinesFilter.Centerlines
            else:
                centerlines = self.Centerlines                
        #write the centerlines if you need to remesh the surface with different parameters (speed up the procedure)
        if self.VTKPrefixFilesName != None:
            if (self.Centerlines == None):
                nomeCenterline = "%s%s" % (self.VTKPrefixFilesName,"_centerlines.vtp")
                self.WritePolyData(centerlines,nomeCenterline)
        self.ComputedCenterline = centerlines

        #FIRST PART FOR BOTH: surface capping and remeshing (only if SkipCapAndRemesh = 0)
        while meshing == "n":
            if self.SurfaceInteractive == 0:
                if self.RemeshingMethod == 'constant':
                    method = "1"
                elif self.RemeshingMethod == 'radius':
                    method = "2"
                elif self.RemeshingMethod == 'edgelengtharray':
                    method = 'edgelengtharray'
                else:
                    self.PrintError('Error: I need constant or radius for remeshing method ')
            elif self.SurfaceInteractive == 1:
                print "\nChoose the interface meshing method:\n"
                method = raw_input("1 = constant size, 2 = radius dependent size  ")
            if method == "1":
                #Set remeshing parameter
                if self.SurfaceInteractive == 0:
                    TargetEdgeLength = self.TargetEdgeLength
                    print "remeshing with constant edgelength: ", TargetEdgeLength
                elif self.SurfaceInteractive == 1:
                    TargetEdgeLength = input("set the edgelength parameter: ")
                #Capping surface
                print "Capping surface"
                capper = vmtkscripts.vmtkSurfaceCapper()
                capper.Surface = surface
                capper.CellEntityIdsArrayName = CellEntityIdsArrayName
                capper.Interactive = 0
                capper.Method = 'simple'
                capper.TriangleOutput = 0
                capper.CellEntityIdOffset = 1
                capper.Execute()
                #Remeshing
                if self.RemeshCapsOnly == 1:
                    remeshing = vmtkscripts.vmtkSurfaceRemeshing()
                    remeshing.Surface = capper.Surface
                    remeshing.CellEntityIdsArrayName = capper.CellEntityIdsArrayName
                    remeshing.TargetEdgeLength = TargetEdgeLength
                    remeshing.ElementSizeMode = 'edgelength'
                    remeshing.ExcludeEntityIds = [1]
                    remeshing.PreserveBoundaryEdges = 1
                    remeshing.Execute()
                else :
                    remeshing = vmtkscripts.vmtkSurfaceRemeshing()
                    remeshing.Surface = capper.Surface
                    remeshing.CellEntityIdsArrayName = capper.CellEntityIdsArrayName
                    remeshing.TargetEdgeLength = TargetEdgeLength
                    remeshing.ElementSizeMode = 'edgelength'
                    remeshing.Execute()
                projection = vmtkscripts.vmtkSurfaceProjection()
                projection.Surface = remeshing.Surface
                projection.ReferenceSurface = capper.Surface
                projection.Execute()
                print "Number of point: ", remeshing.Surface.GetNumberOfPoints()
                print "Number of cells: ", remeshing.Surface.GetNumberOfCells()
                if self.SurfaceInteractive == 0:
                    meshing = "y"
                elif self.SurfaceInteractive == 1:
                    self.OutputText('Displaying.\n')
                    surfaceViewer = vmtkscripts.vmtkSurfaceViewer()
                    surfaceViewer.Surface = remeshing.Surface
                    surfaceViewer.Execute()
                    meshing = raw_input("\nAccept result? (y (yes), n (no)) ")
            elif method == "2":
                #Set remeshing parameter
                if self.SurfaceInteractive == 0:
                    alpha = self.TargetAlphaFactor
                    beta = self.TargetBetaFactor
                    print "remeshing with h =", alpha, " * radius^", beta
                elif self.SurfaceInteractive == 1:
                    print "remeshing with h = alpha * radius^beta"
                    alpha = input("set the alpha factor (multiplicative): ")
                    beta = input("set the beta factor (power): ")
                surfaceDistance = self.ComputeSmartRadius(surface,centerlines,CellEntityIdsArrayName,RadiusArrayName,DistanceToCenterlinesArrayName,SmartDistanceArrayName)
                print "Capping surface"
                capper = vmtkscripts.vmtkSurfaceCapper()
                capper.Surface = surfaceDistance
                capper.CellEntityIdsArrayName=CellEntityIdsArrayName
                capper.Interactive = 0
                capper.Method = 'simple'
                capper.TriangleOutput = 0
                capper.CellEntityIdOffset = 1
                capper.Execute()
                distToCenterlines = capper.Surface.GetPointData().GetArray(SmartDistanceArrayName)
                edgeLengthArray = vtk.vtkFloatArray()
                edgeLengthArray.SetName('EdgeLengthFunction')
                edgeLengthArray.SetNumberOfComponents(1)
                edgeLengthArray.SetNumberOfTuples(capper.Surface.GetNumberOfPoints())
                capper.Surface.GetPointData().AddArray(edgeLengthArray)    
                #funzione target per h = alpha * R^beta
                for j in range(capper.Surface.GetNumberOfPoints()):
                    edgeLengthArray.SetTuple1(j,alpha*math.pow(distToCenterlines.GetComponent(j,0),beta))
                #Remeshing
                if self.RemeshCapsOnly == 1:
                    remeshing = vmtkscripts.vmtkSurfaceRemeshing()
                    remeshing.Surface = capper.Surface
                    remeshing.CellEntityIdsArrayName = capper.CellEntityIdsArrayName
                    remeshing.TargetEdgeLengthArrayName = 'EdgeLengthFunction'
                    remeshing.ElementSizeMode = 'edgelengtharray'
                    remeshing.ExcludeEntityIds = [1]
                    remeshing.PreserveBoundaryEdges = 1
                    remeshing.Execute()
                else :
                    remeshing = vmtkscripts.vmtkSurfaceRemeshing()
                    remeshing.Surface = capper.Surface
                    remeshing.CellEntityIdsArrayName = capper.CellEntityIdsArrayName
                    remeshing.TargetEdgeLengthArrayName = 'EdgeLengthFunction'
                    remeshing.ElementSizeMode = 'edgelengtharray'
                    remeshing.Execute()
                projection = vmtkscripts.vmtkSurfaceProjection()
                projection.Surface = remeshing.Surface
                projection.ReferenceSurface = capper.Surface
                projection.Execute()
                print "Number of point: ", remeshing.Surface.GetNumberOfPoints()
                print "Number of cells: ", remeshing.Surface.GetNumberOfCells()
                if self.SurfaceInteractive == 0:
                    meshing = "y"
                elif self.SurfaceInteractive == 1:
                    self.OutputText('Displaying.\n')
                    surfaceViewer = vmtkscripts.vmtkSurfaceViewer()
                    surfaceViewer.Surface = remeshing.Surface
                    surfaceViewer.Execute()
                    meshing = raw_input("\nAccept result? (y (yes), n (no)) ")
            elif method == 'edgelengtharray':
                if (surface.GetPointData().GetArray(self.TargetEdgeLengthArrayName) == None):
                    self.PrintError('Error: no edgelength array with name specified')
                print "Capping surface"
                capper = vmtkscripts.vmtkSurfaceCapper()
                capper.Surface = surface
                capper.CellEntityIdsArrayName=CellEntityIdsArrayName
                capper.Interactive = 0
                capper.Method = 'simple'
                capper.TriangleOutput = 0
                capper.CellEntityIdOffset = 1
                capper.Execute()
                #Remeshing
                if self.RemeshCapsOnly == 1:
                    remeshing = vmtkscripts.vmtkSurfaceRemeshing()
                    remeshing.Surface = capper.Surface
                    remeshing.CellEntityIdsArrayName = capper.CellEntityIdsArrayName
                    remeshing.TargetEdgeLengthArrayName = self.TargetEdgeLengthArrayName
                    remeshing.ElementSizeMode = 'edgelengtharray'
                    remeshing.ExcludeEntityIds = [1]
                    remeshing.PreserveBoundaryEdges = 1
                    remeshing.Execute()
                else:
                    remeshing = vmtkscripts.vmtkSurfaceRemeshing()
                    remeshing.Surface = capper.Surface
                    remeshing.CellEntityIdsArrayName = capper.CellEntityIdsArrayName
                    remeshing.TargetEdgeLengthArrayName = self.TargetEdgeLengthArrayName
                    remeshing.ElementSizeMode = 'edgelengtharray'
                    remeshing.Execute()
                projection = vmtkscripts.vmtkSurfaceProjection()
                projection.Surface = remeshing.Surface
                projection.ReferenceSurface = capper.Surface
                projection.Execute()
                print "Number of point: ", remeshing.Surface.GetNumberOfPoints()
                print "Number of cells: ", remeshing.Surface.GetNumberOfCells()
                meshing = "y"
            else:
                print "Wrong selected method!"
                meshing = "n"
        
        if  (self.SkipCapAndRemesh == 0):
            #writes the surface with the labels if you need to re-generate the mesh of volume or wall keeping the same interface (speed up the procedure)
            if self.VTKPrefixFilesName != None:
                nomeSurfLabel = "%s%s" % (self.VTKPrefixFilesName,"_labeledsurface.vtp")
                self.WritePolyData(remeshing.Surface,nomeSurfLabel)
            surfaceRemeshed = projection.Surface
            ILabel = 1
        
        self.RemeshedSurface = surfaceRemeshed
        #normals on the whole surface mesh
        normalsFilter = vmtkscripts.vmtkSurfaceNormals()
        normalsFilter.Surface = surfaceRemeshed
        normalsFilter.NormalsArrayName = 'Normals'
        normalsFilter.Execute() 
        #The interface surface is closed by caps. Labels are 1 on the interface and 2,3, ... on the caps.
        #Interface surface extraction:   
        print "interface mesh extraction"
        threshold = vtk.vtkThreshold()
        threshold.SetInputData(surfaceRemeshed)
        threshold.ThresholdBetween(ILabel-0.5,ILabel+0.5)
        threshold.SetInputArrayToProcess(0, 0, 0, 1, CellEntityIdsArrayName)
        threshold.Update()  
        meshToSurface = vmtkscripts.vmtkMeshToSurface()
        meshToSurface.Mesh = threshold.GetOutput()
        meshToSurface.Execute()
        #normals on the interface surface mesh
        inormalsFilter = vmtkscripts.vmtkSurfaceNormals()
        inormalsFilter.Surface = meshToSurface.Surface
        inormalsFilter.NormalsArrayName = 'Normals'
        inormalsFilter.Execute()
        interfaceSurface = inormalsFilter.Surface
        #set interface cell id = InterfaceId
        interfaceCellEntityIdsArray = vtk.vtkIntArray()
        interfaceCellEntityIdsArray.SetName(CellEntityIdsArrayName)
        interfaceCellEntityIdsArray.SetNumberOfTuples(interfaceSurface.GetNumberOfCells())
        interfaceCellEntityIdsArray.FillComponent(0,InterfaceId)
        interfaceSurface.GetCellData().AddArray(interfaceCellEntityIdsArray)
        #set interface point id = InterfaceId
        interfacePointEntityIdsArray = vtk.vtkIntArray()
        interfacePointEntityIdsArray.SetName("PointEntityIdsArray")
        interfacePointEntityIdsArray.SetNumberOfTuples(interfaceSurface.GetNumberOfPoints())
        interfacePointEntityIdsArray.FillComponent(0,InterfaceId)
        interfaceSurface.GetPointData().AddArray(interfacePointEntityIdsArray)
        #boundary interface extraction
        print "boundary extraction"
        boundaryExtractor = vtkvmtk.vtkvmtkPolyDataBoundaryExtractor()
        boundaryExtractor.SetInputData(interfaceSurface)
        boundaryExtractor.Update()
        #-----------recompute normals on rings--------------
        normalCArray = vtk.vtkFloatArray()
        normalCArray.SetName('normalscorrected')
        normalCArray.SetNumberOfComponents(3)
        normalCArray.SetNumberOfTuples(interfaceSurface.GetNumberOfPoints())
        normalCArray.DeepCopy(interfaceSurface.GetPointData().GetNormals())
        interfaceSurface.GetPointData().AddArray(normalCArray)
        point_locator = vtk.vtkPointLocator()
        point_locator.SetDataSet(surfaceRemeshed)
        point_locator.AutomaticOn()
        point_locator.BuildLocator()
        boundaries = vtk.vtkPolyData()
        boundaries = boundaryExtractor.GetOutput()
        pointN = [0.0,0.0,0.0]
        pointNProj = [0.0,0.0,0.0]
        point = [0.0,0.0,0.0]
        pointA = [0.0,0.0,0.0]
        pointB = [0.0,0.0,0.0]
        pointC = [0.0,0.0,0.0]
        planeNormal = [0.0,0.0,0.0]
        n = [0.0,0.0,0.0]
        surfaceRemeshed.GetCellData().SetActiveScalars(CellEntityIdsArrayName)
        areamin = 1000
        areamax = 0
        edgemax = 0
        for i in range (boundaries.GetNumberOfCells()):
            boundary = vtk.vtkPolyLine()
            boundary = boundaries.GetCell(i)
            idb = int(boundaries.GetPointData().GetScalars().GetComponent(boundary.GetPointId(0),0))
            interfaceSurface.GetPoint(idb,pointB)
            idoncomp = point_locator.FindClosestPoint(pointB)
            pointCells = vtk.vtkIdList()
            surfaceRemeshed.GetPointCells(idoncomp,pointCells)
            found = 0
            for k in range (pointCells.GetNumberOfIds()):
                CellInit = pointCells.GetId(k)
                if surfaceRemeshed.GetCellData().GetArray(CellEntityIdsArrayName).GetValue(CellInit) != ILabel and found == 0:
                    labelCap = surfaceRemeshed.GetCellData().GetArray(CellEntityIdsArrayName).GetValue(CellInit)
                    print labelCap
                    cellPoints = vtk.vtkIdList()
                    surfaceRemeshed.GetCellPoints(CellInit,cellPoints)
                    #define the plane
                    surfaceRemeshed.GetPoint(cellPoints.GetId(0),pointA)
                    surfaceRemeshed.GetPoint(cellPoints.GetId(1),pointB)
                    surfaceRemeshed.GetPoint(cellPoints.GetId(2),pointC)
                    triangle = vtk.vtkTriangle()
                    triangle.ComputeNormal(pointA,pointB,pointC,planeNormal)
                    vtk.vtkMath.Normalize(planeNormal)
                    plane = vtk.vtkPlane()
                    plane.SetNormal(planeNormal)
                    plane.SetOrigin(pointB)
                    found = 1
            for j in range (boundary.GetNumberOfPoints()):
                idb = int(boundaries.GetPointData().GetScalars().GetComponent(boundary.GetPointId(j),0))
                interfaceSurface.GetPoint(idb,point)
                pointN = [point[0]+interfaceSurface.GetPointData().GetNormals().GetComponent(idb,0),point[1]+interfaceSurface.GetPointData().GetNormals().GetComponent(idb,1),point[2]+interfaceSurface.GetPointData().GetNormals().GetComponent(idb,2)]
                plane.ProjectPoint(pointN,pointNProj)                 
                n[0] = pointNProj[0] - point[0]
                n[1] = pointNProj[1] - point[1]
                n[2] = pointNProj[2] - point[2]
                vtk.vtkMath.Normalize(n)
                normalCArray.SetTuple3(idb,n[0],n[1],n[2])
            #area and number of cells per cap (this is for the estimate of edgelength on caps)
            threshold = vtk.vtkThreshold()
            threshold.SetInputData(surfaceRemeshed)
            threshold.ThresholdBetween(labelCap-0.5,labelCap+0.5)
            threshold.SetInputArrayToProcess(0, 0, 0, 1, CellEntityIdsArrayName)
            threshold.Update()          
            meshToSurface = vmtkscripts.vmtkMeshToSurface()
            meshToSurface.Mesh = threshold.GetOutput()
            meshToSurface.Execute()
            mass = vtk.vtkMassProperties()
            mass.SetInputData(meshToSurface.Surface)
            mass.Update()
            areaCap = mass.GetSurfaceArea()
            numTri = meshToSurface.Surface.GetNumberOfCells()
            area = areaCap/numTri
            if area < areamin:
                areamin = area
            if area > areamax:
                areamax = area
        edgemax = ((areamax * 4 ) / (3.0**0.5))**0.5
        edgemin = ((areamin * 4 ) / (3.0**0.5))**0.5
        edgemean = edgemax + edgemin / 2

        if self.VTKPrefixFilesName != None:
            nomeOutINTERFACE="%s%s" % (self.VTKPrefixFilesName,"_interface.vtp")
            self.WritePolyData(interfaceSurface,nomeOutINTERFACE)
#-------------------------------------------------------------------------------------------------------------

        if self.Fluid == 1:
            #FLUID MESH GENERATION
            print "FLUID MESH GENERATION"
            if self.FluidInteractive == 0:
                VolumeElementScaleFactor = self.VolumeElementScaleFactor
                if self.BoundaryLayer == 1:
                    bound = "y"
                else:
                    bound = "n"
            elif self.FluidInteractive == 1:
                isok = "n"
                print "volumeelementfactor setted to 0.8"
                isok = raw_input("Type y if you want to change the volumeelementfactor parameter. ")
                if isok == "y":
                    VolumeElementScaleFactor = input("Input the desired volumeelementfactor (e.g. 0.8) " )
                else:
                    VolumeElementScaleFactor = 0.8
                print "Do you want a fluid boundary layer?"
                bound = raw_input("y (yes), n (no). (Default no.) ")
                if bound != "y" and bound != "n":
                    print "No boundary layer"
                    bound = "n"
            if bound == "y":
                placeholderCellEntityId = 9999
                if self.FluidInteractive == 0:
                    NumberOfSubLayers = self.BoundaryLayerNumberOfSubLayers
                    SubLayerRatio = self.SubLayerRatio
                    if self.BoundaryLayerMethod == "constant":
                        bl_method = "1"
                    elif self.BoundaryLayerMethod == "radius":
                        bl_method = "2"
                    else:
                        self.PrintError('Error: I need constant or radius for remeshing method ')
                elif self.FluidInteractive == 1:
                    NumberOfSubLayers = input("Number of desired sublayers: ")
                    SubLayerRatio = input("SubLayerRatio i.e. ratio between the thickness of two successive boundary layers (e.g. 0.8):")
                    print "Choose the desired boundary layer:"
                    bl_method = raw_input("1 = constant boundary layer, 2 = radius dependent boundary layer  ")
                if bl_method == "1":
                    if self.FluidInteractive == 0:
                        BoundaryLayerThickness = self.BoundaryLayerConstantThickness
                    elif self.FluidInteractive == 1:
                        BoundaryLayerThickness = input ("desired boundary layer constant thickness:")
                    surfaceToMesh = vmtkscripts.vmtkSurfaceToMesh()
                    surfaceToMesh.Surface = interfaceSurface
                    surfaceToMesh.Execute()
                    #BoundaryLayer creation
                    boundaryLayer = vmtkscripts.vmtkBoundaryLayer()
                    boundaryLayer.Mesh = surfaceToMesh.Mesh
                    if self.NormalsCorrected == 1:
                        boundaryLayer.WarpVectorsArrayName = 'normalscorrected'
                    else:
                        boundaryLayer.WarpVectorsArrayName = 'Normals'
                    boundaryLayer.NegateWarpVectors = True
                    boundaryLayer.ConstantThickness = True
                    boundaryLayer.NumberOfSubLayers = NumberOfSubLayers
                    boundaryLayer.NumberOfSubsteps = self.BoundaryLayerNumberOfSubsteps
                    boundaryLayer.IncludeSurfaceCells = 0
                    boundaryLayer.IncludeSidewallCells = 0
                    boundaryLayer.SubLayerRatio = SubLayerRatio
                    boundaryLayer.Thickness = BoundaryLayerThickness
                    boundaryLayer.ThicknessRatio = 1.0
                    boundaryLayer.SidewallCellEntityId = placeholderCellEntityId                    
                    boundaryLayer.InnerSurfaceCellEntityId = 1
                    boundaryLayer.Execute()
                elif bl_method == "2":
                    if self.FluidInteractive == 0:
                        alpha2 = self.BoundaryLayerAlphaFactor
                        beta2 = self.BoundaryLayerBetaFactor
                        print "Boundary layer generation with thickness =", alpha2 ," * radius^", beta2
                    elif self.FluidInteractive == 1:
                        print "Boundary layer generation with thickness = alpha2 * radius^beta2"
                        alpha2 = input("set the alpha2 factor (multiplicative): ")
                        beta2 = input("set the beta2 factor (power): ")
                    CDistanceSurface = self.ComputeSmartRadius(interfaceSurface,centerlines,CellEntityIdsArrayName,RadiusArrayName,DistanceToCenterlinesArrayName,SmartDistanceArrayName)
                    thicknessNew = vtk.vtkFloatArray()
                    thicknessNew.SetName('thickness')
                    thicknessNew.SetNumberOfComponents(1)
                    thicknessNew.SetNumberOfTuples(CDistanceSurface.GetNumberOfPoints())
                    CDistanceSurface.GetPointData().AddArray(thicknessNew)
                    DistToCenterlines = CDistanceSurface.GetPointData().GetArray(SmartDistanceArrayName)
                    for j in range(CDistanceSurface.GetNumberOfPoints()):
                        thicknessNew.SetTuple1(j,alpha2*math.pow(DistToCenterlines.GetComponent(j,0),beta2))        
                    #boundary layer generation with thicknessNew    
                    surfaceToMesh = vmtkscripts.vmtkSurfaceToMesh()
                    surfaceToMesh.Surface = CDistanceSurface
                    surfaceToMesh.Execute()
                    boundaryLayer = vmtkscripts.vmtkBoundaryLayer()
                    boundaryLayer.Mesh = surfaceToMesh.Mesh
                    if self.NormalsCorrected == 1:
                        boundaryLayer.WarpVectorsArrayName = 'normalscorrected'
                        print "I am here"
                    else:
                        boundaryLayer.WarpVectorsArrayName = 'Normals'
                    boundaryLayer.NegateWarpVectors = True
                    boundaryLayer.ThicknessArrayName = 'thickness'
                    boundaryLayer.NumberOfSubLayers = NumberOfSubLayers
                    boundaryLayer.NumberOfSubsteps = self.BoundaryLayerNumberOfSubsteps
                    boundaryLayer.IncludeSurfaceCells = 0
                    boundaryLayer.IncludeSidewallCells = 0
                    boundaryLayer.SubLayerRatio = SubLayerRatio
                    boundaryLayer.ThicknessRatio = 1.0
                    boundaryLayer.SidewallCellEntityId = placeholderCellEntityId                    
                    boundaryLayer.InnerSurfaceCellEntityId = 1 #wallEntityOffset
                    boundaryLayer.Execute()  
                #Fill boundaryLayer cellEntityIdsArray     
                cellEntityIdsArray = vtk.vtkIntArray()
                cellEntityIdsArray.SetName(CellEntityIdsArrayName)
                cellEntityIdsArray.SetNumberOfTuples(boundaryLayer.Mesh.GetNumberOfCells())
                cellEntityIdsArray.FillComponent(0,0.0)
                boundaryLayer.Mesh.GetCellData().AddArray(cellEntityIdsArray)
                innerCellEntityIdsArray = vtk.vtkIntArray()
                innerCellEntityIdsArray.SetName(CellEntityIdsArrayName)
                innerCellEntityIdsArray.SetNumberOfTuples(boundaryLayer.InnerSurfaceMesh.GetNumberOfCells())
                innerCellEntityIdsArray.FillComponent(0,400) #the inner surface is just support. impose a high value id.
                boundaryLayer.InnerSurfaceMesh.GetCellData().AddArray(innerCellEntityIdsArray)
                meshToSurface = vmtkscripts.vmtkMeshToSurface()
                meshToSurface.Mesh = boundaryLayer.InnerSurfaceMesh
                meshToSurface.Execute()           
                innerSurface = meshToSurface.Surface                     
                #inner capping:
                print "Capping inner"
                capper = vmtkscripts.vmtkSurfaceCapper()
                capper.Surface = innerSurface
                capper.CellEntityIdsArrayName=CellEntityIdsArrayName
                capper.Interactive = 0
                capper.Method = 'simple'
                capper.TriangleOutput = 0
                capper.CellEntityIdOffset = 1
                capper.Execute()
                print "remeshing caps"
                if self.SkipCapAndRemesh == 1:
                    if self.FluidInteractive == 0:
                        if self.CapRemeshingMethod == "constant":
                            cap_method = "1"
                            CapEdgeLength = self.CapTargetEdgeLength
                        elif self.CapRemeshingMethod == "radius":
                            cap_method = "2"
                            alpha3 = self.CapTargetAlphaFactor
                            beta3 = self.CapTargetBetaFactor
                            print "Cap generation with h = ", alpha3 ," * radius ^ ", beta3
                    elif self.FluidInteractive == 1:
                        print "Choose the method to remesh caps:"
                        cap_method = raw_input("1 = constant cap remeshing, 2 = radius dependent cap remeshing ")
                    if cap_method == "1":
                        if self.FluidInteractive == 1:
                            print "the mean edge length of your original input caps is", edgemean, "(max edgelength=",edgemax, "min edgelength=", edgemin,")"
                            CapEdgeLength = input ("desired edgelength for caps:")
                        remeshing = vmtkscripts.vmtkSurfaceRemeshing()
                        remeshing.Surface = capper.Surface
                        remeshing.CellEntityIdsArrayName = capper.CellEntityIdsArrayName
                        remeshing.TargetEdgeLength = CapEdgeLength
                        remeshing.ElementSizeMode = 'edgelength'
                        #remeshing.MaxEdgeLength = edgemax
                        remeshing.ExcludeEntityIds = [400]
                        remeshing.PreserveBoundaryEdges = 1
                        remeshing.Execute()
                    else:
                        if self.FluidInteractive == 1:
                            print "Cap generation with h = alpha * radius^beta"
                            alpha3 = input("set the alpha factor (multiplicative): ")
                            beta3 = input("set the beta factor (power): ")                    
                        capSurface = self.ComputeSmartRadius(capper.Surface,centerlines,CellEntityIdsArrayName,RadiusArrayName,DistanceToCenterlinesArrayName,SmartDistanceArrayName)
                        CapEdgeLengthArray = vtk.vtkFloatArray()
                        CapEdgeLengthArray.SetName('CapEdgeLengthFunction')
                        CapEdgeLengthArray.SetNumberOfComponents(1)
                        CapEdgeLengthArray.SetNumberOfTuples(capper.Surface.GetNumberOfPoints())
                        capSurface.GetPointData().AddArray(CapEdgeLengthArray)
                        DistToCenterlines = capSurface.GetPointData().GetArray(SmartDistanceArrayName)
                        #funzione target per h = alpha * R^beta
                        for j in range(capSurface.GetNumberOfPoints()):
                            CapEdgeLengthArray.SetTuple1(j, alpha3 * math.pow(DistToCenterlines.GetComponent(j,0),beta3))
                        #Remeshing
                        remeshing = vmtkscripts.vmtkSurfaceRemeshing()
                        remeshing.Surface = capSurface
                        remeshing.CellEntityIdsArrayName = CellEntityIdsArrayName
                        remeshing.TargetEdgeLengthArrayName = 'CapEdgeLengthFunction'
                        remeshing.ElementSizeMode = 'edgelengtharray'
                        remeshing.ExcludeEntityIds = [400]
                        remeshing.PreserveBoundaryEdges = 1
                        remeshing.Execute()                                       
                else:
                    print ""
                    if method == "1":
                        print "I remesh caps with constant edgelength equal to", TargetEdgeLength
                        remeshing = vmtkscripts.vmtkSurfaceRemeshing()
                        remeshing.Surface = capper.Surface
                        remeshing.CellEntityIdsArrayName = capper.CellEntityIdsArrayName
                        remeshing.TargetEdgeLength = TargetEdgeLength
                        remeshing.ElementSizeMode = 'edgelength'
                        remeshing.MaxEdgeLength = edgemax      
                        remeshing.ExcludeEntityIds = [400]
                        remeshing.PreserveBoundaryEdges = 1
                        remeshing.Execute()
                    elif method == "2":
                        print "I remesh caps with h =", alpha ,"* radius^", beta
                        capSurface = self.ComputeSmartRadius(capper.Surface,centerlines,CellEntityIdsArrayName,RadiusArrayName,DistanceToCenterlinesArrayName,SmartDistanceArrayName)
                        CapEdgeLengthArray = vtk.vtkFloatArray()
                        CapEdgeLengthArray.SetName('CapEdgeLengthFunction')
                        CapEdgeLengthArray.SetNumberOfComponents(1)
                        CapEdgeLengthArray.SetNumberOfTuples(capper.Surface.GetNumberOfPoints())
                        capSurface.GetPointData().AddArray(CapEdgeLengthArray)
                        DistToCenterlines = capSurface.GetPointData().GetArray(SmartDistanceArrayName)
                        #funzione target per h = alpha * R^beta
                        for j in range(capSurface.GetNumberOfPoints()):
                            CapEdgeLengthArray.SetTuple1(j, alpha * math.pow(DistToCenterlines.GetComponent(j,0),beta))
                        #Remeshing
                        remeshing = vmtkscripts.vmtkSurfaceRemeshing()
                        remeshing.Surface = capSurface
                        remeshing.CellEntityIdsArrayName = CellEntityIdsArrayName
                        remeshing.TargetEdgeLengthArrayName = 'CapEdgeLengthFunction'
                        remeshing.ElementSizeMode = 'edgelengtharray'
                        remeshing.ExcludeEntityIds = [400]
                        remeshing.PreserveBoundaryEdges = 1
                        remeshing.Execute()
                    elif method == "edgelengtharray":
                        projection = vmtkscripts.vmtkSurfaceProjection()
                        projection.Surface = capper.Surface
                        projection.ReferenceSurface = surfaceRemeshed
                        projection.Execute()
                        #Remeshing
                        remeshing = vmtkscripts.vmtkSurfaceRemeshing()
                        remeshing.Surface = capper.Surface
                        remeshing.CellEntityIdsArrayName = CellEntityIdsArrayName
                        remeshing.TargetEdgeLengthArrayName = self.TargetEdgeLengthArrayName
                        remeshing.ElementSizeMode = 'edgelengtharray'
                        remeshing.ExcludeEntityIds = [400]
                        remeshing.PreserveBoundaryEdges = 1
                        remeshing.Execute()
                #Append interface + mesh boundary layer
                appendFilterBL = vtkvmtk.vtkvmtkAppendFilter() 
                appendFilterBL.AddInputData(interfaceSurface)
                appendFilterBL.AddInputData(boundaryLayer.Mesh)
                appendFilterBL.Update()            
                #tetra the result
                tetrahedralizeBL = vtk.vtkDataSetTriangleFilter()
                tetrahedralizeBL.SetInputData(appendFilterBL.GetOutput())
                tetrahedralizeBL.Update()
                #this is only for region growing!
                appendFilterBL2 = vtkvmtk.vtkvmtkAppendFilter() 
                appendFilterBL2.AddInputData(innerSurface)
                appendFilterBL2.AddInputData(tetrahedralizeBL.GetOutput())
                appendFilterBL2.Update()          
                #define the BL surface
                meshToSurface = vmtkscripts.vmtkMeshToSurface()
                meshToSurface.Mesh = appendFilterBL2.GetOutput()
                meshToSurface.Execute()
                surfaceBL = meshToSurface.Surface         
                #threshold to obtain a mesh of only inner caps
                threshold = vtk.vtkThreshold()
                threshold.SetInputData(remeshing.Surface)
                threshold.ThresholdByLower(100)
                threshold.SetInputArrayToProcess(0, 0, 0, 1, CellEntityIdsArrayName)
                threshold.Update()
                meshToSurface = vmtkscripts.vmtkMeshToSurface()
                meshToSurface.Mesh = threshold.GetOutput()
                meshToSurface.Execute()
                surfaceCaps = vtk.vtkPolyData()
                surfaceCaps = meshToSurface.Surface            
                #locator on Caps to find correct id for rings
                point_locatorCaps = vtk.vtkPointLocator()
                point_locatorCaps.SetDataSet(surfaceCaps)
                point_locatorCaps.AutomaticOn()
                point_locatorCaps.BuildLocator()
                #CellIdArray on Caps
                capsCellIdArray = surfaceCaps.GetCellData().GetArray(CellEntityIdsArrayName)        
                #locator on surface BL for the region growing on caps
                point_locatorBL = vtk.vtkPointLocator()
                point_locatorBL.SetDataSet(surfaceBL)
                point_locatorBL.AutomaticOn()
                point_locatorBL.BuildLocator()
                #CellIdArray on surface BL          
                BLCellEntityIdsArray = surfaceBL.GetCellData().GetArray(CellEntityIdsArrayName)  
                #boundary extractor on inner surface (boundaries for the region growing on rings)
                BLboundaryExtractor = vtkvmtk.vtkvmtkPolyDataBoundaryExtractor()
                BLboundaryExtractor.SetInputData(innerSurface)
                BLboundaryExtractor.Update()
                print "start region growing on internal rings"
                FlagIdsArray = vtk.vtkIntArray()
                FlagIdsArray.SetName("Flag")
                FlagIdsArray.SetNumberOfTuples(surfaceBL.GetNumberOfCells())
                FlagIdsArray.FillComponent(0,0.0)
                surfaceBL.GetCellData().AddArray(FlagIdsArray)
                pointB = [0.0,0.0,0.0]
                boundaries = vtk.vtkPolyData()
                boundaries = BLboundaryExtractor.GetOutput()
                numberOfEdges = 0
                for i in range (boundaries.GetNumberOfCells()):
                    boundary = vtk.vtkPolyLine()
                    boundary = boundaries.GetCell(i)
                    idb = int(boundaries.GetPointData().GetScalars().GetComponent(boundary.GetPointId(1),0))
                    innerSurface.GetPoint(idb,pointB)
                    idoncaps = point_locatorCaps.FindClosestPoint(pointB)  
                    pointCellsCaps = vtk.vtkIdList()
                    surfaceCaps.GetPointCells(idoncaps,pointCellsCaps)
                    RingId = capsCellIdArray.GetComponent(pointCellsCaps.GetId(0),0)           
                    print "RingId equal to", RingId
                    idvero = point_locatorBL.FindClosestPoint(pointB)  
                    pointCells = vtk.vtkIdList()
                    surfaceBL.GetPointCells(idvero,pointCells)
                    found = 0
                    #choice of the first cell
                    for k in range (pointCells.GetNumberOfIds()):
                        CellInit = pointCells.GetId(k)
                        if BLCellEntityIdsArray.GetValue(CellInit) == 0 and found == 0:
                            cellPoints = vtk.vtkIdList()
                            surfaceBL.GetCellPoints(CellInit,cellPoints)
                            for j in range (0,3):
                                if cellPoints.GetId(j) != idvero:
                                    cellPointsList = vtk.vtkIdList()
                                    cellPointsList.InsertNextId(idvero)
                                    cellPointsList.InsertNextId(cellPoints.GetId(j))
                                    cellNeighborList = vtk.vtkIdList()
                                    surfaceBL.GetCellNeighbors(CellInit,cellPointsList,cellNeighborList)
                                    if BLCellEntityIdsArray.GetValue(cellNeighborList.GetId(0)) != 0:
                                        pointRot = idvero
                                        pointApp = cellPoints.GetId(j)
                                        for p in range(0,3):
                                            if cellPoints.GetId(p) != pointRot and cellPoints.GetId(p) != pointApp:
                                                pointMove = cellPoints.GetId(p)
                                        Cell = CellInit
                                        found = 1
                    if found == 0:
                        print "error! no element boundary found"
                    else:
                        #print "RingId", RingId
                        BLCellEntityIdsArray.SetTuple1(Cell,RingId)
                        FlagIdsArray.SetTuple1(Cell,1)
                    check = 1
                    layer = 1
                    while check == 1:
                        RotEdge = vtk.vtkIdList()
                        RotEdge.InsertNextId(pointRot)
                        RotEdge.InsertNextId(pointMove)
                        cellNeighborList = vtk.vtkIdList()
                        surfaceBL.GetCellNeighbors(Cell,RotEdge,cellNeighborList)
                        CellNeighbor = cellNeighborList.GetId(0)
                        if BLCellEntityIdsArray.GetValue(CellNeighbor) == 0:
                            BLCellEntityIdsArray.SetTuple1(CellNeighbor,RingId)
                            FlagIdsArray.SetTuple1(CellNeighbor,layer)
                            cellPoints = vtk.vtkIdList()
                            surfaceBL.GetCellPoints(CellNeighbor,cellPoints)
                            for j in range (0,3):
                                iduse = cellPoints.GetId(j)
                                if iduse != pointRot and iduse != pointMove:
                                    pointMove2 = iduse
                            pointMove = pointMove2
                            Cell = CellNeighbor
                        else:
                            if (BLCellEntityIdsArray.GetValue(CellNeighbor) != 0) and (BLCellEntityIdsArray.GetValue(CellNeighbor) != RingId):
                                cellPoints = vtk.vtkIdList()
                                surfaceBL.GetCellPoints(Cell,cellPoints)
                                for j in range (0,3):
                                    iduse = cellPoints.GetId(j)
                                    if iduse != pointRot and iduse != pointMove:
                                        pointMove2 = iduse
                                pointRot = pointMove
                                pointMove = pointMove2
                                Cell = Cell
                            elif BLCellEntityIdsArray.GetValue(CellNeighbor) == RingId and FlagIdsArray.GetValue(CellNeighbor) == layer - 1:
                                cellPoints = vtk.vtkIdList()
                                surfaceBL.GetCellPoints(Cell,cellPoints)
                                for j in range (0,3):
                                    iduse = cellPoints.GetId(j)
                                    if iduse != pointRot and iduse != pointMove:
                                        pointMove2 = iduse
                                pointRot = pointMove
                                pointMove = pointMove2
                                Cell = Cell
                            elif BLCellEntityIdsArray.GetValue(CellNeighbor) == RingId and FlagIdsArray.GetValue(CellNeighbor) == layer:
                                pointCells = vtk.vtkIdList()
                                surfaceBL.GetPointCells(pointMove,pointCells)
                                found = 0
                                for k in range (pointCells.GetNumberOfIds()):
                                    CellInit = pointCells.GetId(k)
                                    if BLCellEntityIdsArray.GetValue(CellInit) == 0 and found == 0:
                                        cellPoints = vtk.vtkIdList()
                                        surfaceBL.GetCellPoints(CellInit,cellPoints)
                                        for j in range (0,3):
                                            iduse = cellPoints.GetId(j)
                                            if cellPoints.GetId(j) != pointMove:
                                                cellPointsList = vtk.vtkIdList()
                                                cellPointsList.InsertNextId(pointMove)
                                                cellPointsList.InsertNextId(cellPoints.GetId(j))
                                                cellNeighborList = vtk.vtkIdList()
                                                surfaceBL.GetCellNeighbors(CellInit,cellPointsList,cellNeighborList)
                                                if BLCellEntityIdsArray.GetValue(cellNeighborList.GetId(0))==RingId:
                                                    pointRot = pointMove
                                                    pointApp = cellPoints.GetId(j)
                                                    for p in range(0,3):
                                                        if cellPoints.GetId(p) != pointRot and cellPoints.GetId(p) != pointApp:
                                                            pointMove = cellPoints.GetId(p)
                                                    Cell = CellInit
                                                    found = 1     
                                if found == 0:
                                    check = 0
                                else:
                                    BLCellEntityIdsArray.SetTuple1(Cell,RingId)
                                    layer = layer+1
                                    FlagIdsArray.SetTuple1(Cell,layer)                                  
                #threshold surface BL to erease the inner surface
                threshold = vtk.vtkThreshold()
                threshold.SetInputData(surfaceBL)
                threshold.ThresholdByLower(290)
                threshold.SetInputArrayToProcess(0, 0, 0, 1, CellEntityIdsArrayName)
                threshold.Update()
                meshToSurface = vmtkscripts.vmtkMeshToSurface()
                meshToSurface.Mesh = threshold.GetOutput()
                meshToSurface.Execute()
                surfaceBLext = meshToSurface.Surface                        
                #fluid volume mesh geneartion inside the boundaryLayer (tetgen)
                print "Generating fluid volume mesh with TetGen:"
                sizingFunction = vtkvmtk.vtkvmtkPolyDataSizingFunction()
                sizingFunction.SetInputData(remeshing.Surface)
                sizingFunction.SetSizingFunctionArrayName(SizingFunctionArrayName)
                sizingFunction.SetScaleFactor(VolumeElementScaleFactor)
                sizingFunction.Update()
                surfaceToMesh = vmtkscripts.vmtkSurfaceToMesh()
                surfaceToMesh.Surface = sizingFunction.GetOutput()
                surfaceToMesh.Execute()
                #tetgen meshing         
                tetgen = vmtkscripts.vmtkTetGen()
                tetgen.Mesh = surfaceToMesh.Mesh
                tetgen.GenerateCaps = 0
                tetgen.UseSizingFunction = 1
                tetgen.SizingFunctionArrayName = SizingFunctionArrayName
                tetgen.CellEntityIdsArrayName = CellEntityIdsArrayName
                tetgen.Order = 1
                tetgen.Quality = 1
                tetgen.PLC = 1
                tetgen.NoBoundarySplit = 1
                tetgen.RemoveSliver = 1
                tetgen.OutputSurfaceElements = 0
                tetgen.OutputVolumeElements = 1 
                tetgen.Execute()#extract the "tetra part" from tetrahedralizeBL and append with the "surface part"
                fluidTetra = tetrahedralizeBL.GetOutput()
                cta = vtk.vtkIntArray()
                cta.SetName("celltypes")
                cta.SetNumberOfComponents(1)
                cta.SetNumberOfTuples(fluidTetra.GetNumberOfCells())
                fluidTetra.GetCellData().AddArray(cta)
                for i in range(fluidTetra.GetNumberOfCells()):
                    cta.SetTuple1(i,fluidTetra.GetCell(i).GetCellType())
                thre = vtk.vtkThreshold()
                thre.SetInputData(fluidTetra)
                thre.ThresholdBetween(10-1,10+1)
                thre.SetInputArrayToProcess(0, 0, 0, 1, "celltypes")
                thre.Update()               
                appendFilter = vtkvmtk.vtkvmtkAppendFilter()
                appendFilter.AddInputData(surfaceBLext)
                appendFilter.AddInputData(thre.GetOutput())
                appendFilter.AddInputData(tetgen.Mesh)
                appendFilter.AddInputData(surfaceCaps)
                appendFilter.Update()
                #final tetrahedralization (results = volumetric mesh with all interfaces)
                tetrahedralizeFluid = vtk.vtkDataSetTriangleFilter()
                tetrahedralizeFluid.SetInputData(appendFilter.GetOutput())
                tetrahedralizeFluid.Update()
                Fluid = tetrahedralizeFluid.GetOutput()
                print "fluid cells", Fluid.GetNumberOfCells()
                print "fluid points", Fluid.GetNumberOfPoints()
            elif bound == "n":
                sizingFunction = vtkvmtk.vtkvmtkPolyDataSizingFunction()
                sizingFunction.SetInputData(normalsFilter.Surface)
                sizingFunction.SetSizingFunctionArrayName(SizingFunctionArrayName)
                sizingFunction.SetScaleFactor(VolumeElementScaleFactor)
                sizingFunction.Update()           
                surfaceToMesh2 = vmtkscripts.vmtkSurfaceToMesh()
                surfaceToMesh2.Surface = sizingFunction.GetOutput()
                surfaceToMesh2.Execute()
                #tetgen meshing
                tetgen = vmtkscripts.vmtkTetGen()
                tetgen.Mesh = surfaceToMesh2.Mesh
                tetgen.GenerateCaps = 0
                tetgen.UseSizingFunction = 1
                tetgen.SizingFunctionArrayName = SizingFunctionArrayName
                tetgen.CellEntityIdsArrayName = CellEntityIdsArrayName
                tetgen.Order = 1
                tetgen.Quality = 1
                tetgen.PLC = 1
                tetgen.NoBoundarySplit = 1
                tetgen.RemoveSliver = 1
                tetgen.OutputSurfaceElements = 1
                tetgen.OutputVolumeElements = 1
                tetgen.Execute()
                #final tetrahedralization (results = volumetric mesh with all interfaces)
                tetrahedralizeFluid = vtk.vtkDataSetTriangleFilter()
                tetrahedralizeFluid.SetInputData(tetgen.Mesh)
                tetrahedralizeFluid.Update()
                #Inteface surface mesh id change (note, by default is equal to one, we set it equal to InterfaceId)
                Fluid = tetrahedralizeFluid.GetOutput()
                FluidcellIdsArray = Fluid.GetCellData().GetArray(CellEntityIdsArrayName)
                for j in range(Fluid.GetNumberOfCells()):
                    if FluidcellIdsArray.GetComponent(j,0) == 1:
                        FluidcellIdsArray.SetTuple1(j,InterfaceId)
            print "edge definition"
            boundaries = vtk.vtkPolyData()
            boundaries = boundaryExtractor.GetOutput()
            numberOfEdges = 0
            for i in range (boundaries.GetNumberOfCells()):
                boundary = vtk.vtkPolyLine()
                boundary = boundaries.GetCell(i)
                numberOfEdges = numberOfEdges + boundary.GetNumberOfPoints()
            EdgeCellArray1 = vtk.vtkIntArray()
            EdgeCellArray1.SetName("EdgeCellArray1")
            EdgeCellArray1.SetNumberOfTuples(numberOfEdges)
            EdgeCellArray1.SetNumberOfComponents(3)
            EdgeCellArray1.FillComponent(0,0)
            EdgeCellArray1.FillComponent(1,0)
            EdgeCellArray1.FillComponent(2,0)
            Fluid.BuildLinks()
            cellEntityIdsArray = vtk.vtkIntArray()
            cellEntityIdsArray.DeepCopy(Fluid.GetCellData().GetArray(CellEntityIdsArrayName))
            pointEntityIdsArray = vtk.vtkIntArray()
            pointEntityIdsArray.SetName("WritingPointEntityIdsArray")
            pointEntityIdsArray.SetNumberOfTuples(Fluid.GetNumberOfPoints())
            pointEntityIdsArray.FillComponent(0,0)
            Fluid.GetPointData().AddArray(pointEntityIdsArray)
            for i in range(Fluid.GetNumberOfPoints()):
                point = Fluid.GetPoint(i)
                pointCells = vtk.vtkIdList()
                Fluid.GetPointCells(i,pointCells)
                minTriangleCellEntityId = -1
                tetraCellEntityId = -1
                for j in range(pointCells.GetNumberOfIds()):
                    cellId = pointCells.GetId(j)
                    if Fluid.GetCellType(cellId) == triangleCellType:
                        cellEntityId = cellEntityIdsArray.GetValue(cellId)
                        if cellEntityId < minTriangleCellEntityId or minTriangleCellEntityId == -1:
                            minTriangleCellEntityId = cellEntityId
                    else:
                        tetraCellEntityId = cellEntityIdsArray.GetValue(cellId)
                cellEntityId = tetraCellEntityId
                if minTriangleCellEntityId != -1:
                    cellEntityId = minTriangleCellEntityId
                pointEntityIdsArray.SetTuple1(i,cellEntityId)
            point_locator = vtk.vtkPointLocator()
            point_locator.SetDataSet(Fluid)
            point_locator.AutomaticOn()
            point_locator.BuildLocator()
            pointB = [0.0,0.0,0.0]
            index = 0
            for i in range (boundaries.GetNumberOfCells()):
                boundary = vtk.vtkPolyLine()
                boundary = boundaries.GetCell(i)
                idb = int(boundaries.GetPointData().GetScalars().GetComponent(boundary.GetPointId(1),0))
                interfaceSurface.GetPoint(idb,pointB)
                idvero = point_locator.FindClosestPoint(pointB)  
                pointCells = vtk.vtkIdList()
                Fluid.GetPointCells(idvero,pointCells)
                for j in range(pointCells.GetNumberOfIds()):
                    cellId = pointCells.GetId(j)
                    if Fluid.GetCellType(cellId) == triangleCellType:
                        cellEntityId = cellEntityIdsArray.GetValue(cellId)
                        if cellEntityId != InterfaceId:
                            EdgeCellEntityId = cellEntityId * 10
                boundary = boundaries.GetCell(i)
                for j in range (boundary.GetNumberOfPoints()):
                    idb = int(boundaries.GetPointData().GetScalars().GetComponent(boundary.GetPointId(j),0))
                    interfaceSurface.GetPoint(idb,pointB)     
                    idbvero = point_locator.FindClosestPoint(pointB)
                    interfacePointEntityIdsArray.SetTuple1(idb,EdgeCellEntityId)
                    pointEntityIdsArray.SetTuple1(idbvero,EdgeCellEntityId)
                    if j == (boundary.GetNumberOfPoints() - 1):
                        idb2 = int(boundaries.GetPointData().GetScalars().GetComponent(boundary.GetPointId(0),0))
                        interfaceSurface.GetPoint(idb2,pointB)
                        idb2vero = point_locator.FindClosestPoint(pointB)
                        EdgeCellArray1.InsertTuple3(index,idbvero + 1, idb2vero + 1, EdgeCellEntityId)
                        index = index + 1
                    else:
                        idb2 = int(boundaries.GetPointData().GetScalars().GetComponent(boundary.GetPointId(j+1),0))
                        interfaceSurface.GetPoint(idb2,pointB)
                        idb2vero = point_locator.FindClosestPoint(pointB)
                        EdgeCellArray1.InsertTuple3(index,idbvero+1, idb2vero+1, EdgeCellEntityId)
                        index = index + 1
            if self.VTKPrefixFilesName != None:
                nomeOutMeshVTU = "%s%s" % (self.VTKPrefixFilesName,"_FluidMesh.vtu")
                self.WriteUnstructuredGrid(Fluid,nomeOutMeshVTU)
            self.FluidMesh = Fluid
            #Write LifeV file with edges
            print "writing LifeV mesh"
            OutputFileName = self.LifeVFluidMeshName
            f = open(OutputFileName, 'w')
            line = "MeshVersionFormatted 1\n\n"
            line += "Dimension\n"
            line += "3\n\n"
            line += "Vertices\n"
            line += "%d\n" % Fluid.GetNumberOfPoints()
            f.write(line)
            for i in range(Fluid.GetNumberOfPoints()):
                point = Fluid.GetPoint(i)
                pointEntityId = pointEntityIdsArray.GetValue(i)
                line = "%f  %f  %f  %d\n" % (point[0]*self.ScaleFactor, point[1]*self.ScaleFactor, point[2]*self.ScaleFactor, pointEntityId)
                f.write(line)
            line = "\n"
            #tetras
            tetraCellIdArray = vtk.vtkIdTypeArray()
            Fluid.GetIdsOfCellsOfType(tetraCellType,tetraCellIdArray)
            numberOfTetras = tetraCellIdArray.GetNumberOfTuples()
            line += "Tetrahedra\n"
            line += "%d\n" % numberOfTetras
            f.write(line)
            for i in range(numberOfTetras):
                tetraCellId = tetraCellIdArray.GetValue(i) 
                cellPointIds = Fluid.GetCell(tetraCellId).GetPointIds()
                line = ''
                for j in range(cellPointIds.GetNumberOfIds()):
                    if j > 0:
                        line += '  '
                    line += "%d" % (cellPointIds.GetId(j)+1)
                cellEntityId = cellEntityIdsArray.GetValue(tetraCellId)   
                line += '  %d\n' % (cellEntityId+1)
                f.write(line)
            line = "\n"
            #tirangles
            triangleCellIdArray = vtk.vtkIdTypeArray()
            Fluid.GetIdsOfCellsOfType(triangleCellType,triangleCellIdArray)
            numberOfTriangles = triangleCellIdArray.GetNumberOfTuples()
            line += "Triangles\n"
            line += "%d\n" % numberOfTriangles
            f.write(line)
            for i in range(numberOfTriangles):
                triangleCellId = triangleCellIdArray.GetValue(i)
                cellPointIds = Fluid.GetCell(triangleCellId).GetPointIds()
                line = ''
                for j in range(cellPointIds.GetNumberOfIds()):
                    if j > 0:
                        line += '  '
                    line += "%d" % (cellPointIds.GetId(j)+1)
                cellEntityId = cellEntityIdsArray.GetValue(triangleCellId)
                line += '  %d\n' % cellEntityId
                f.write(line)
            line = "\n"
            #Edges  
            line += "Edges\n"
            line += "%d\n" % numberOfEdges
            f.write(line)
            for i in range (numberOfEdges):
                line = "%d  %d  %d\n" % (EdgeCellArray1.GetComponent(i,0), EdgeCellArray1.GetComponent(i,1), EdgeCellArray1.GetComponent(i,2))
                f.write(line)
            f.close()
            print 'Fluid LifeV mesh file written: ', self.LifeVFluidMeshName
            print ' '

        elif self.Fluid == 0:
            #edge definition for solid construction 
            print "edge definition"
            boundaries = vtk.vtkPolyData()
            boundaries = boundaryExtractor.GetOutput()
            cellEntityIdsArray = vtk.vtkIntArray()
            cellEntityIdsArray.DeepCopy(surfaceRemeshed.GetCellData().GetArray(CellEntityIdsArrayName))
            point_locator = vtk.vtkPointLocator()
            point_locator.SetDataSet(surfaceRemeshed)
            point_locator.AutomaticOn()
            point_locator.BuildLocator()
            pointB = [0.0,0.0,0.0]
            for i in range (boundaries.GetNumberOfCells()):
                boundary = vtk.vtkPolyLine()
                boundary = boundaries.GetCell(i)
                idb = int(boundaries.GetPointData().GetScalars().GetComponent(boundary.GetPointId(1),0))
                interfaceSurface.GetPoint(idb,pointB)
                idvero = point_locator.FindClosestPoint(pointB)  
                pointCells = vtk.vtkIdList()
                surfaceRemeshed.GetPointCells(idvero,pointCells)
                for j in range(pointCells.GetNumberOfIds()):
                    cellId = pointCells.GetId(j)
                    cellEntityId = cellEntityIdsArray.GetValue(cellId)
                    if cellEntityId != ILabel:
                        EdgeCellEntityId = cellEntityId * 10
                boundary = boundaries.GetCell(i)
                for j in range (boundary.GetNumberOfPoints()):
                    idb = int(boundaries.GetPointData().GetScalars().GetComponent(boundary.GetPointId(j),0))
                    interfacePointEntityIdsArray.SetTuple1(idb,EdgeCellEntityId)

##-------------------------------------------------------------------------------------------
        #SOLID MESH GENERATION
        if self.Solid:
            print "STRUCTURE MESH GENERATION"
            #---parameters initialization----
            WallVolumePtsId = 0.0

            surfaceDistance = self.ComputeSmartRadius(interfaceSurface,centerlines,CellEntityIdsArrayName,RadiusArrayName,DistanceToCenterlinesArrayName,SmartDistanceArrayName)
            surfMesh = vmtkscripts.vmtkSurfaceToMesh()
            surfMesh.Surface = surfaceDistance
            surfMesh.Execute()
            InterfaceMesh = vtk.vtkUnstructuredGrid()
            InterfaceMesh = surfMesh.Mesh
            #thickness array
            thickness = vtk.vtkFloatArray()
            thickness.SetName('thicknessFactor')
            thickness.SetNumberOfComponents(1)
            thickness.SetNumberOfTuples(InterfaceMesh.GetNumberOfPoints())
            InterfaceMesh.GetPointData().AddArray(thickness)
            #id arrays
            InterfacePointEntityIdsArray = vtk.vtkIntArray()
            InterfacePointEntityIdsArray.SetName('InterfacePointEntityIds')
            InterfacePointEntityIdsArray.SetNumberOfTuples(InterfaceMesh.GetNumberOfPoints())
            InterfacePointEntityIdsArray.FillComponent(0,InterfaceId)
            InterfaceMesh.GetPointData().AddArray(InterfacePointEntityIdsArray)
            InterfaceCellEntityIdsArray = vtk.vtkIntArray()
            InterfaceCellEntityIdsArray.SetName(CellEntityIdsArrayName)
            InterfaceCellEntityIdsArray.SetNumberOfTuples(InterfaceMesh.GetNumberOfCells())
            InterfaceCellEntityIdsArray.FillComponent(0,InterfaceId)
            InterfaceMesh.GetCellData().AddArray(InterfaceCellEntityIdsArray)
            #extrusion parameters
            cellPts = vtk.vtkIdList()
            DistCenterline = InterfaceMesh.GetPointData().GetArray(SmartDistanceArrayName)
            if self.SolidInteractive == 0:
                SolidNumberOfSubLayers = self.SolidNumberOfSubLayers
                if self.SolidMethod == "radius":
                    extrusion_method = "2"
                    thicknessRatioC = self.SolidAlphaFactor
                elif self.SolidMethod == "constant":
                    extrusion_method = "1"
                    constant_thickness = self.SolidConstantThickness
            elif self.SolidInteractive == 1:
                SolidNumberOfSubLayers = input ("Number of sublayers desired to represent the structure: ")
                print "Choose the extrusion method for the principal vessel:"
                extrusion_method = raw_input("1 = constant extrusion, 2 = radius dependent extrusion  ")
            if extrusion_method == "1":
                if self.SolidInteractive == 1:
                    constant_thickness = input ("input the desired constant thickness: ")
                for j in range(InterfaceMesh.GetNumberOfPoints()):
                    thickness.SetTuple1(j,constant_thickness)
            elif extrusion_method == "2":
                if self.SolidInteractive == 1:
                    thicknessRatioC = input ("input the desired alpha factor ( thickess = alpha * r):  (e.g. 0.1) ")
                for j in range(InterfaceMesh.GetNumberOfPoints()):
                    thickness.SetTuple1(j, DistCenterline.GetComponent(j,0) * thicknessRatioC)  
            else:
                print "Attention, option", extrusion_method, " does not exist! I'll extrude the area with a radius dependent extrusion method and a thickness ratio of 0.2"
                thicknessRatioC = 0.2
                for j in range(InterfaceMesh.GetNumberOfPoints()):
                    thickness.SetTuple1(j, DistCenterline.GetComponent(j,0) * thicknessRatioC)

            MultiArea = "n"
            if self.SolidInteractive == 1:
                MultiArea = raw_input("Do you want to delineate multiple areas for different extrusion parameters and/or volume conditions? (y (yes), n (no)). Default is no. ")
                #multi area step. only for different volumes and/or different extrusions
            if MultiArea == "y":
                region_delineation = 1
                InnerRegionId = 0
            else:
                region_delineation = 0
                InnerRegionId = 0
            while region_delineation == 1:
                InnerRegionId = InnerRegionId + 1
                InnerRegionIdToSet = InnerRegionId + 100
                print "Please delineate the area", InnerRegionId, "on the mesh (left click + ctrl ) "
                Meshsurf = vmtkscripts.vmtkMeshToSurface()
                InterfaceMesh.GetCellData().Update()
                InterfaceMesh.GetPointData().Update()
                Meshsurf.Mesh = InterfaceMesh
                Meshsurf.Execute()
                Meshsurf.Surface.GetPointData().Update()
                Meshsurf.Surface.GetPointData().SetActiveScalars('InterfacePointEntityIds')
                lut = vtk.vtkLookupTable()
                lut.SetTableRange (110,200)
                mapper = vtk.vtkPolyDataMapper() 
                mapper.SetInputData(Meshsurf.Surface)
                mapper.ScalarVisibilityOn()
                mapper.SetScalarModeToUsePointData()
                mapper.SetColorModeToMapScalars()
                mapper.UseLookupTableScalarRangeOn() 
                mapper.SetLookupTable(lut)
                mapper.Update()
                actor = vtk.vtkActor()
                actor.SetMapper(mapper)
                actor.GetProperty().SetInterpolationToFlat()
                actor.GetProperty().SetRepresentationToSurface()
                actor.GetProperty().SetEdgeVisibility(1)
                actor.GetProperty().SetEdgeColor(0 , 0, 0)
                renderer = vtk.vtkRenderer()
                renderer.AddActor(actor)
                renderer.SetBackground(1, 1, 1)
                renderWindow = vtk.vtkRenderWindow()
                renderWindow.AddRenderer(renderer)
                interactor = vtk.vtkRenderWindowInteractor()
                interactorStyle = vtk.vtkInteractorStyleTrackballCamera()
                interactor.SetInteractorStyle(interactorStyle)
                interactor.SetRenderWindow(renderWindow)
                contourWidget = vtk.vtkContourWidget()
                contourWidget.SetInteractor(interactor)
                rep = vtk.vtkOrientedGlyphContourRepresentation()
                rep = contourWidget.GetRepresentation()
                rep.GetLinesProperty().SetColor(1, 0, 0)
                rep.GetLinesProperty().SetLineWidth(10.0)
                pointPlacer = vtk.vtkPolygonalSurfacePointPlacer()
                pointPlacer.AddProp(actor)
                pointPlacer.GetPolys().AddItem (Meshsurf.Surface)
                rep.SetPointPlacer(pointPlacer)
                polycontour = vtk.vtkPolyData()
                polycontour = rep.GetContourRepresentationAsPolyData()
                interpolator = vtk.vtkPolygonalSurfaceContourLineInterpolator()
                interpolator.GetPolys().AddItem(Meshsurf.Surface)
                rep.SetLineInterpolator(interpolator) 
                renderWindow.Render()
                interactor.Initialize()
                contourWidget.EnabledOn()
                interactor.Start() 
                select = vtk.vtkSelectPolyData()
                select.SetInputData(Meshsurf.Surface)
                select.SetLoop(polycontour.GetPoints())
                region = raw_input("Do you want the Smallest (1) or The Biggest (2) region? (Default Smallest)")
                if region == "2":
                    select.SetSelectionModeToLargestRegion()
                else:
                    select.SetSelectionModeToSmallestRegion()
                select.GenerateSelectionScalarsOn()
                select.Update()
                InsideOutsideRegion = select.GetOutput().GetPointData().GetScalars()  
                print "Choose the extrusion method for the selected vessel area:"
                extrusion_method = raw_input("1 = constant extrusion, 2 = radius dependent extrusion  ")
                count = 0
                if extrusion_method == "1":
                    constant_thickness = input ("input desired constant thickness: ")
                    for j in range(InterfaceMesh.GetNumberOfPoints()):
                        if InsideOutsideRegion.GetComponent(j,0) <= 0:
                            thickness.SetTuple1(j,constant_thickness)
                            InterfacePointEntityIdsArray.SetTuple1(j,InnerRegionIdToSet)
                            count = count + 1
                elif extrusion_method == "2":
                    thicknessRatioC = input ("input desired alpha factor ( thickess = alpha * r):  (e.g. 0.1) ")
                    for j in range(InterfaceMesh.GetNumberOfPoints()):
                        if InsideOutsideRegion.GetComponent(j,0) <= 0:
                            thickness.SetTuple1(j,DistCenterline.GetComponent(j,0)*thicknessRatioC)  
                            InterfacePointEntityIdsArray.SetTuple1(j,InnerRegionIdToSet)
                            count = count + 1
                else:
                    print "Attention, option", extrusion_method, " does not exist! I'll treat the selected area as the principal area"   
                count = 0
                for j in range (InterfaceMesh.GetNumberOfCells()):
                    InterfaceMesh.GetCellPoints(j,cellPts)
                    if(InterfacePointEntityIdsArray.GetComponent(cellPts.GetId(0),0) == InnerRegionIdToSet) and  (InterfacePointEntityIdsArray.GetComponent(cellPts.GetId(1),0) == InnerRegionIdToSet) and (InterfacePointEntityIdsArray.GetComponent(cellPts.GetId(2),0) == InnerRegionIdToSet):
                        InterfaceCellEntityIdsArray.SetTuple1(j,InnerRegionIdToSet)
                        count = count + 1
                print "Number of selected cells:", count
                region_delineation = raw_input("Do you want to delineate another region? (y (yes), n (no)) ")
                if region_delineation == "y":
                    region_delineation = 1
                else:
                    region_delineation = 0

            print "Surface extrusion"
            InterfaceMesh = self.SmoothPointDataSimple(InterfaceMesh,'thicknessFactor',1)
            Extrusion = vmtkscripts.vmtkBoundaryLayer()
            Extrusion.Mesh = InterfaceMesh
            if self.NormalsCorrected == 1:
                Extrusion.WarpVectorsArrayName = 'normalscorrected'
            else:
                Extrusion.WarpVectorsArrayName = 'Normals'
            Extrusion.ThicknessArrayName = 'thicknessFactor'
            Extrusion.NumberOfSubLayers = SolidNumberOfSubLayers
            Extrusion.ThicknessRatio = 1.0
            Extrusion.ConstantThickness = 0
            Extrusion.IncludeSurfaceCells = 0
            Extrusion.IncludeSidewallCells = 0
            Extrusion.NumberOfSubsteps = self.SolidNumberOfSubsteps
            Extrusion.Execute()
            ExtrusioncellEntityIdsArray = vtk.vtkIntArray()
            ExtrusioncellEntityIdsArray.SetName(CellEntityIdsArrayName)
            ExtrusioncellEntityIdsArray.SetNumberOfTuples(Extrusion.Mesh.GetNumberOfCells())
            ExtrusioncellEntityIdsArray.FillComponent(0,WallVolumePtsId)
            Extrusion.Mesh.GetCellData().AddArray(ExtrusioncellEntityIdsArray)
            meshToSurfaceExtrusionInner = vmtkscripts.vmtkMeshToSurface()
            meshToSurfaceExtrusionInner.Mesh = Extrusion.InnerSurfaceMesh
            meshToSurfaceExtrusionInner.Execute()
            surfaceExternal = vtk.vtkPolyData()
            surfaceExternal = meshToSurfaceExtrusionInner.Surface 
            #cell array on surface external
            OuterPointEntityIdsArray = vtk.vtkIntArray()
            OuterPointEntityIdsArray.SetName(CellEntityIdsArrayName)
            OuterPointEntityIdsArray.SetNumberOfTuples(surfaceExternal.GetNumberOfPoints())
            OuterPointEntityIdsArray.FillComponent(0,ExternalId)
            surfaceExternal.GetPointData().AddArray(OuterPointEntityIdsArray)
            OuterCellEntityIdsArray = vtk.vtkIntArray()
            OuterCellEntityIdsArray.SetName(CellEntityIdsArrayName)
            OuterCellEntityIdsArray.SetNumberOfTuples(surfaceExternal.GetNumberOfCells())
            OuterCellEntityIdsArray.FillComponent(0,ExternalId)
            surfaceExternal.GetCellData().AddArray(OuterCellEntityIdsArray)
            #multiple region definitions. only for different boundary conditions
            MultiAreaB = "n"
            if self.SolidInteractive == 1:
                MultiAreaB = raw_input("Do you want to delineate multiple areas for different bundary conditions (areas not equal to part A delineation)? (y (yes), n (no)). Default is no. ")
            if MultiAreaB == "y":
                region_delineation = 1
            else:
                region_delineation = 0
            while region_delineation == 1:
                InnerRegionId = InnerRegionId + 1
                OuterRegionId = ExternalId + InnerRegionId * 10
                print "Please delineate the area", InnerRegionId, "on the mesh (left click + ctrl ) "
                surfaceExternal.GetCellData().Update()
                surfaceExternal.GetPointData().Update()
                surfaceExternal.GetPointData().SetActiveScalars(CellEntityIdsArrayName)
                lut = vtk.vtkLookupTable()
                lut.SetTableRange (210,250)
                mapper = vtk.vtkPolyDataMapper() 
                mapper.SetInputData(surfaceExternal)
                mapper.ScalarVisibilityOn()
                mapper.SetScalarModeToUsePointData()
                mapper.SetColorModeToMapScalars()
                mapper.UseLookupTableScalarRangeOn() 
                mapper.SetLookupTable(lut)
                mapper.Update()
                actor = vtk.vtkActor()
                actor.SetMapper(mapper)
                actor.GetProperty().SetInterpolationToFlat()
                actor.GetProperty().SetRepresentationToSurface()
                actor.GetProperty().SetEdgeVisibility(1)
                actor.GetProperty().SetEdgeColor(0 , 0, 0)
                renderer = vtk.vtkRenderer()
                renderer.AddActor(actor)
                renderer.SetBackground(1, 1, 1)
                renderWindow = vtk.vtkRenderWindow()
                renderWindow.AddRenderer(renderer)
                interactor = vtk.vtkRenderWindowInteractor()
                interactorStyle = vtk.vtkInteractorStyleTrackballCamera()
                interactor.SetInteractorStyle(interactorStyle)
                interactor.SetRenderWindow(renderWindow)
                contourWidget = vtk.vtkContourWidget()
                contourWidget.SetInteractor(interactor)
                rep = vtk.vtkOrientedGlyphContourRepresentation()
                rep = contourWidget.GetRepresentation()
                rep.GetLinesProperty().SetColor(0, 0, 1)
                rep.GetLinesProperty().SetLineWidth(10.0) 
                pointPlacer = vtk.vtkPolygonalSurfacePointPlacer()
                pointPlacer.AddProp(actor)
                pointPlacer.GetPolys().AddItem (surfaceExternal)
                rep.SetPointPlacer(pointPlacer)
                polycontour = vtk.vtkPolyData()
                polycontour = rep.GetContourRepresentationAsPolyData()
                interpolator = vtk.vtkPolygonalSurfaceContourLineInterpolator()
                interpolator.GetPolys().AddItem(surfaceExternal)
                rep.SetLineInterpolator(interpolator) 
                renderWindow.Render()
                interactor.Initialize()
                contourWidget.EnabledOn()
                interactor.Start() 
                select = vtk.vtkSelectPolyData()
                select.SetInputData(surfaceExternal)
                select.SetLoop(polycontour.GetPoints())
                region = raw_input("Do you want the Smallest (1) or The Biggest (2) region? (Default Smallest)")
                if region == "2":
                    select.SetSelectionModeToLargestRegion()
                else:
                    select.SetSelectionModeToSmallestRegion()
                select.GenerateSelectionScalarsOn()
                select.Update()
                InsideOutsideRegion = select.GetOutput().GetPointData().GetScalars()  
                count = 0
                for j in range(surfaceExternal.GetNumberOfPoints()):
                    if InsideOutsideRegion.GetComponent(j,0) <= 0:
                        OuterPointEntityIdsArray.SetTuple1(j,OuterRegionId)
                        count = count + 1
                for j in range (surfaceExternal.GetNumberOfCells()):
                    surfaceExternal.GetCellPoints(j,cellPts)
                    if(OuterPointEntityIdsArray.GetComponent(cellPts.GetId(0),0) == OuterRegionId) and  (OuterPointEntityIdsArray.GetComponent(cellPts.GetId(1),0) == OuterRegionId) and (OuterPointEntityIdsArray.GetComponent(cellPts.GetId(2),0) == OuterRegionId):
                        OuterCellEntityIdsArray.SetTuple1(j,OuterRegionId)                                  
                region_delineation = raw_input("Do you want to delineate another region? (y (yes), n (no)) ")
                if region_delineation == "y":
                    region_delineation = 1
                else:
                    region_delineation = 0
            #Append
            appendFilterSolid = vtkvmtk.vtkvmtkAppendFilter() 
            appendFilterSolid.AddInputData(Extrusion.Mesh)
            appendFilterSolid.AddInputData(InterfaceMesh)
            appendFilterSolid.AddInputData(surfaceExternal)
            appendFilterSolid.Update()
            Solid = appendFilterSolid.GetOutput()
            appendCellIdArray = Solid.GetCellData().GetArray(CellEntityIdsArrayName)
      
            RGCellEntityIdsArray = vtk.vtkIntArray()
            RGCellEntityIdsArray.SetName("RGCellEntityIds")
            RGCellEntityIdsArray.SetNumberOfTuples(Solid.GetNumberOfCells())
            RGCellEntityIdsArray.DeepCopy(appendCellIdArray)
            Solid.GetCellData().AddArray(RGCellEntityIdsArray)

            if (MultiArea == "y" or MultiAreaB == "yes"):
                print "start region growing in volumes"
                triangleCellIdArray = vtk.vtkIdTypeArray()
                Solid.GetIdsOfCellsOfType(triangleCellType,triangleCellIdArray)
                numberOfTriangles = triangleCellIdArray.GetNumberOfTuples()
                pointListTriangle = vtk.vtkIdList()
                pointListTriangleNotEqual = vtk.vtkIdList()
                cellNeighborList = vtk.vtkIdList()
                pointListPenta = vtk.vtkIdList()
                for i in range(numberOfTriangles): #go only on tringles
                    triangleCellId = triangleCellIdArray.GetValue(i)
                    SurfEntityId = appendCellIdArray.GetValue(triangleCellId) #cellId
                    if SurfEntityId > 100 and SurfEntityId < InterfaceId: #search for new ragions
                        AreaId = SurfEntityId-100
                        AreaVolumePtsId = AreaId
                        AreaExternalId = ExternalId+AreaId*10
                        stop = 0
                        Cell = triangleCellId
                        Solid.GetCellPoints(Cell,pointListTriangle)
                        pointListTriangleNotEqual.DeepCopy(pointListTriangle)
                        Solid.GetCellNeighbors(Cell,pointListTriangleNotEqual,cellNeighborList)                    
                        if cellNeighborList.GetNumberOfIds() != 1:
                            print "too much cells"
                        while stop == 0:
                            Solid.GetCellNeighbors(Cell,pointListTriangleNotEqual,cellNeighborList)                    
                            if cellNeighborList.GetNumberOfIds() != 1:
                                print "there's a problem"
                            if appendCellIdArray.GetValue(cellNeighborList.GetId(0)) == ExternalId: #if I am on the external face
                                appendCellIdArray.SetTuple1(cellNeighborList.GetId(0),AreaExternalId)
                                appendCellIdArray.SetTuple1(triangleCellId,InterfaceId) #reassign the value on the interface
                                stop = 1
                            elif appendCellIdArray.GetValue(cellNeighborList.GetId(0)) == WallVolumePtsId: #if I am in the volume
                                appendCellIdArray.SetTuple1(cellNeighborList.GetId(0),AreaVolumePtsId)
                                Cell = cellNeighborList.GetId(0)
                                Solid.GetCellPoints(Cell,pointListPenta)
                                pointListTriangleNotEqual.Reset()
                                for j in range(pointListPenta.GetNumberOfIds()):
                                    Id1 = pointListPenta.GetId(j)
                                    equal = 0
                                    for k in range (0,3):
                                        if pointListTriangle.GetId(k) == Id1:
                                            equal = 1
                                    if equal == 0:
                                        pointListTriangleNotEqual.InsertNextId(Id1)
                                pointListTriangle.DeepCopy(pointListTriangleNotEqual)
            #tetra
            tetrahedralize = vtk.vtkDataSetTriangleFilter()
            tetrahedralize.SetInputData(Solid)
            tetrahedralize.Update()
            #surface extraction
            meshToSurfaceExtrusion = vmtkscripts.vmtkMeshToSurface()
            meshToSurfaceExtrusion.Mesh = tetrahedralize.GetOutput()
            meshToSurfaceExtrusion.Execute()
            surfaceTotale = vtk.vtkPolyData()
            surfaceTotale = meshToSurfaceExtrusion.Surface
            #DA VERIFICARE. non ricordo piu perche ho messo due array di questo tipo? mi servono?
            TotalCellEntityIdsArray = surfaceTotale.GetCellData().GetArray(CellEntityIdsArrayName)
            RGCellEntityIdsArray = surfaceTotale.GetCellData().GetArray(CellEntityIdsArrayName)
            #RGCellEntityIdsArray = surfaceTotale.GetCellData().GetArray("RGCellEntityIds")
            print "start region growing on rings"
            #Region growing to assign different ids to rings
            FlagIdsArray = vtk.vtkIntArray()
            FlagIdsArray.SetName("Flag")
            FlagIdsArray.SetNumberOfTuples(surfaceTotale.GetNumberOfCells())
            FlagIdsArray.FillComponent(0, 0.0)
            surfaceTotale.GetCellData().AddArray(FlagIdsArray)
            print surfaceTotale.GetNumberOfCells()
            point_locator = vtk.vtkPointLocator()
            point_locator.SetDataSet(surfaceTotale)
            point_locator.AutomaticOn()
            point_locator.BuildLocator()
            pointB = [0.0,0.0,0.0]
            boundaries = vtk.vtkPolyData()
            boundaries = boundaryExtractor.GetOutput()
            numberOfEdges = 0
            for i in range (boundaries.GetNumberOfCells()):
                boundary = vtk.vtkPolyLine()
                boundary = boundaries.GetCell(i)
                numberOfEdges = numberOfEdges+boundary.GetNumberOfPoints()
            for i in range (boundaries.GetNumberOfCells()):
                boundary = vtk.vtkPolyLine()
                boundary = boundaries.GetCell(i)
                idb = int(boundaries.GetPointData().GetScalars().GetComponent(boundary.GetPointId(1),0))
                RingId = interfacePointEntityIdsArray.GetComponent(idb,0) / 10
                print "RingId equal to", RingId
                interfaceSurface.GetPoint(idb,pointB)
                idvero = point_locator.FindClosestPoint(pointB)  
                pointCells = vtk.vtkIdList()
                surfaceTotale.GetPointCells(idvero,pointCells)
                found = 0
                for k in range (pointCells.GetNumberOfIds()):
                    CellInit = pointCells.GetId(k)
                    if RGCellEntityIdsArray.GetValue(CellInit) == 0 and found == 0:
                        cellPoints = vtk.vtkIdList()
                        surfaceTotale.GetCellPoints(CellInit,cellPoints)
                        for j in range (0,3):
                            if cellPoints.GetId(j) != idvero:
                                cellPointsList = vtk.vtkIdList()
                                cellPointsList.InsertNextId(idvero)
                                cellPointsList.InsertNextId(cellPoints.GetId(j))
                                cellNeighborList = vtk.vtkIdList()
                                surfaceTotale.GetCellNeighbors(CellInit,cellPointsList,cellNeighborList)
                                if RGCellEntityIdsArray.GetValue(cellNeighborList.GetId(0)) != 0:
                                    pointRot = idvero
                                    pointApp = cellPoints.GetId(j)
                                    for p in range(0,3):
                                        if cellPoints.GetId(p) != pointRot and cellPoints.GetId(p) != pointApp:
                                            pointMove = cellPoints.GetId(p)
                                    Cell = CellInit
                                    found = 1
                if found == 0:
                    print "error! no element boundary found"
                else:
                    TotalCellEntityIdsArray.SetTuple1(Cell,RingId)
                    RGCellEntityIdsArray.SetTuple1(Cell,RingId)
                    FlagIdsArray.SetTuple1(Cell,1)
                check = 1
                layer = 1
                while check == 1:
                    RotEdge = vtk.vtkIdList()
                    RotEdge.InsertNextId(pointRot)
                    RotEdge.InsertNextId(pointMove)
                    cellNeighborList = vtk.vtkIdList()
                    surfaceTotale.GetCellNeighbors(Cell,RotEdge,cellNeighborList)
                    CellNeighbor = cellNeighborList.GetId(0)
                    if RGCellEntityIdsArray.GetValue(CellNeighbor) == 0:
                        TotalCellEntityIdsArray.SetTuple1(CellNeighbor,RingId)
                        RGCellEntityIdsArray.SetTuple1(CellNeighbor,RingId)
                        FlagIdsArray.SetTuple1(CellNeighbor,layer)
                        cellPoints = vtk.vtkIdList()
                        surfaceTotale.GetCellPoints(CellNeighbor,cellPoints)
                        for j in range (0,3):
                            iduse = cellPoints.GetId(j)
                            if iduse != pointRot and iduse != pointMove:
                                pointMove2 = iduse
                        pointMove = pointMove2
                        Cell = CellNeighbor
                    else:
                        if (RGCellEntityIdsArray.GetValue(CellNeighbor) != 0) and (RGCellEntityIdsArray.GetValue(CellNeighbor) != RingId):
                            cellPoints = vtk.vtkIdList()
                            surfaceTotale.GetCellPoints(Cell,cellPoints)
                            for j in range (0,3):
                                iduse = cellPoints.GetId(j)
                                if iduse != pointRot and iduse != pointMove:
                                    pointMove2 = iduse
                            pointRot = pointMove
                            pointMove = pointMove2
                            Cell = Cell
                        elif RGCellEntityIdsArray.GetValue(CellNeighbor) == RingId and FlagIdsArray.GetValue(CellNeighbor) == layer - 1:
                            cellPoints = vtk.vtkIdList()
                            surfaceTotale.GetCellPoints(Cell,cellPoints)
                            for j in range (0,3):
                                iduse = cellPoints.GetId(j)
                                if iduse != pointRot and iduse != pointMove:
                                    pointMove2 = iduse
                            pointRot = pointMove
                            pointMove = pointMove2
                            Cell = Cell
                        elif RGCellEntityIdsArray.GetValue(CellNeighbor) == RingId and FlagIdsArray.GetValue(CellNeighbor) == layer:
                            pointCells = vtk.vtkIdList()
                            surfaceTotale.GetPointCells(pointMove,pointCells)
                            found = 0
                            for k in range (pointCells.GetNumberOfIds()):
                                CellInit = pointCells.GetId(k)
                                if RGCellEntityIdsArray.GetValue(CellInit) == 0 and found == 0:
                                    cellPoints = vtk.vtkIdList()
                                    surfaceTotale.GetCellPoints(CellInit,cellPoints)
                                    for j in range (0,3):
                                        iduse = cellPoints.GetId(j)
                                        if cellPoints.GetId(j) != pointMove:
                                            cellPointsList = vtk.vtkIdList()
                                            cellPointsList.InsertNextId(pointMove)        
                                            cellPointsList.InsertNextId(cellPoints.GetId(j))
                                            cellNeighborList = vtk.vtkIdList()
                                            surfaceTotale.GetCellNeighbors(CellInit,cellPointsList,cellNeighborList)
                                            if RGCellEntityIdsArray.GetValue(cellNeighborList.GetId(0)) == RingId:
                                                pointRot = pointMove
                                                pointApp = cellPoints.GetId(j)
                                                for p in range(0,3):
                                                    if cellPoints.GetId(p) != pointRot and cellPoints.GetId(p) != pointApp:
                                                        pointMove = cellPoints.GetId(p)
                                                Cell = CellInit
                                                found = 1     
                            if found == 0:
                                check = 0
                            else:
                                TotalCellEntityIdsArray.SetTuple1(Cell,RingId)
                                RGCellEntityIdsArray.SetTuple1(Cell,RingId)
                                layer = layer + 1
                                FlagIdsArray.SetTuple1(Cell,layer)        
            #create the tetra part and the surface part and append them
            solidTetra = tetrahedralize.GetOutput()
            cta = vtk.vtkIntArray()
            cta.SetName("celltypes")
            cta.SetNumberOfComponents(1)
            cta.SetNumberOfTuples(solidTetra.GetNumberOfCells())
            solidTetra.GetCellData().AddArray(cta)
            for i in range(solidTetra.GetNumberOfCells()):
                cta.SetTuple1(i,solidTetra.GetCell(i).GetCellType())
            thre = vtk.vtkThreshold()
            thre.SetInputData(solidTetra)
            thre.ThresholdBetween(10-1,10+1)
            thre.SetInputArrayToProcess(0, 0, 0, 1, "celltypes")
            thre.Update()       
            #Append surfaceTotale and the volumetric mesh
            appendFilterSolid = vtkvmtk.vtkvmtkAppendFilter() 
            appendFilterSolid.AddInputData(surfaceTotale)
            appendFilterSolid.AddInputData(thre.GetOutput())
            appendFilterSolid.Update()
            #create the final solid
            Solid = appendFilterSolid.GetOutput()
            Solid.BuildLinks()       
    ##-----------------------------------------------------------------------
            #Edge definition: external boundary extraction
            #External surface mesh creation 
            meshToSurfaceExternal = vmtkscripts.vmtkMeshToSurface()
            meshToSurfaceExternal.Mesh = Extrusion.InnerSurfaceMesh
            meshToSurfaceExternal.Execute()
            boundaryExtractorExt = vtkvmtk.vtkvmtkPolyDataBoundaryExtractor()
            boundaryExtractorExt.SetInputData(meshToSurfaceExternal.Surface)
            boundaryExtractorExt.Update()
            boundariesExt = vtk.vtkPolyData()
            boundariesExt = boundaryExtractorExt.GetOutput()
            numberOfEdgesExt = 0
            for i in range (boundariesExt.GetNumberOfCells()):
                boundaryExt = vtk.vtkPolyLine()
                boundaryExt = boundariesExt.GetCell(i)
                numberOfEdgesExt = numberOfEdgesExt + boundaryExt.GetNumberOfPoints()
            numberOfEdgesTot = numberOfEdgesExt + numberOfEdges
            cellEntityIdsArray = vtk.vtkIntArray()
            cellEntityIdsArray.DeepCopy(Solid.GetCellData().GetArray(CellEntityIdsArrayName))
            #final array
            SolidpointEntityIdsArray = vtk.vtkIntArray()
            SolidpointEntityIdsArray.SetName("WritingPointEntityIdsArray")
            SolidpointEntityIdsArray.SetNumberOfTuples(Solid.GetNumberOfPoints())
            SolidpointEntityIdsArray.FillComponent(0,0)
            Solid.GetPointData().AddArray(SolidpointEntityIdsArray)
            #edge array and definition for external edges
            EdgeCellArray2 = vtk.vtkIntArray()
            EdgeCellArray2.SetName("EdgeCellArray1")
            EdgeCellArray2.SetNumberOfTuples(numberOfEdgesExt)
            EdgeCellArray2.SetNumberOfComponents(3)
            EdgeCellArray2.FillComponent(0,0)
            EdgeCellArray2.FillComponent(1,0)
            EdgeCellArray2.FillComponent(2,0)
            #I set the tetraCellEntityId all equal to 1 because LifeV doesn't support different values for volume vertex!
            for i in range(Solid.GetNumberOfPoints()):
                point = Solid.GetPoint(i)
                pointCells = vtk.vtkIdList()
                Solid.GetPointCells(i,pointCells)
                minTriangleCellEntityId = -1
                tetraCellEntityId = -1
                for j in range(pointCells.GetNumberOfIds()):
                    cellId = pointCells.GetId(j)
                    if Solid.GetCellType(cellId) == triangleCellType:
                        cellEntityId = cellEntityIdsArray.GetValue(cellId)
                        if cellEntityId < minTriangleCellEntityId or minTriangleCellEntityId == -1:
                            minTriangleCellEntityId = cellEntityId
                    else:
                        tetraCellEntityId = 0
                cellEntityId = tetraCellEntityId
                if minTriangleCellEntityId != -1:
                    cellEntityId = minTriangleCellEntityId
                SolidpointEntityIdsArray.SetTuple1(i,cellEntityId)
            point_locator = vtk.vtkPointLocator()
            point_locator.SetDataSet(Solid)
            point_locator.AutomaticOn()
            point_locator.BuildLocator()
            pointB = [0.0,0.0,0.0]
            #edge array and definition for internal edges
            EdgeCellArray1 = vtk.vtkIntArray()
            EdgeCellArray1.SetName("EdgeCellArray1")
            EdgeCellArray1.SetNumberOfTuples(numberOfEdges)
            EdgeCellArray1.SetNumberOfComponents(3)
            EdgeCellArray1.FillComponent(0,0)
            EdgeCellArray1.FillComponent(1,0)
            EdgeCellArray1.FillComponent(2,0)
            index = 0
            for i in range (boundaries.GetNumberOfCells()):
                boundary = vtk.vtkPolyLine()
                boundary = boundaries.GetCell(i)
                idb = int(boundaries.GetPointData().GetScalars().GetComponent(boundary.GetPointId(1),0))
                interfaceSurface.GetPoint(idb,pointB)
                idvero = point_locator.FindClosestPoint(pointB)  
                pointCells = vtk.vtkIdList()
                Solid.GetPointCells(idvero,pointCells)
                for j in range(pointCells.GetNumberOfIds()):
                    cellId = pointCells.GetId(j)
                    if Solid.GetCellType(cellId) == triangleCellType:
                        cellEntityId = cellEntityIdsArray.GetValue(cellId)
                        if cellEntityId != InterfaceId:
                            EdgeCellEntityId = cellEntityId * 10
                boundary = boundaries.GetCell(i)
                for j in range (boundary.GetNumberOfPoints()):
                    idb = int(boundaries.GetPointData().GetScalars().GetComponent(boundary.GetPointId(j),0))
                    interfaceSurface.GetPoint(idb,pointB)
                    idbvero = point_locator.FindClosestPoint(pointB)
                    SolidpointEntityIdsArray.SetTuple1(idbvero,EdgeCellEntityId)
                    if j == (boundary.GetNumberOfPoints()-1):
                        idb2 = int(boundaries.GetPointData().GetScalars().GetComponent(boundary.GetPointId(0),0))
                        interfaceSurface.GetPoint(idb2,pointB)
                        idb2vero = point_locator.FindClosestPoint(pointB)
                        EdgeCellArray1.InsertTuple3(index,idbvero+1, idb2vero+1, EdgeCellEntityId)
                        index = index + 1
                    else:
                        idb2 = int(boundaries.GetPointData().GetScalars().GetComponent(boundary.GetPointId(j+1),0))
                        interfaceSurface.GetPoint(idb2,pointB)
                        idb2vero = point_locator.FindClosestPoint(pointB)
                        EdgeCellArray1.InsertTuple3(index,idbvero+1, idb2vero+1, EdgeCellEntityId)
                        index = index + 1
            index = 0
            for i in range (boundariesExt.GetNumberOfCells()):
                boundary = vtk.vtkPolyLine()
                boundary = boundariesExt.GetCell(i)
                idb = int(boundariesExt.GetPointData().GetScalars().GetComponent(boundary.GetPointId(1),0))
                meshToSurfaceExternal.Surface.GetPoint(idb,pointB)
                idvero = point_locator.FindClosestPoint(pointB)  
                pointCells = vtk.vtkIdList()
                Solid.GetPointCells(idvero,pointCells)
                for j in range(pointCells.GetNumberOfIds()):
                    cellId = pointCells.GetId(j)
                    if Solid.GetCellType(cellId) == triangleCellType:
                        cellEntityId = cellEntityIdsArray.GetValue(cellId)
                        if cellEntityId < 200:
                            EdgeCellEntityId = (cellEntityId * 10) + 1
                boundary = boundariesExt.GetCell(i)
                for j in range (boundary.GetNumberOfPoints()):    
                    idb = int(boundariesExt.GetPointData().GetScalars().GetComponent(boundary.GetPointId(j),0))
                    meshToSurfaceExternal.Surface.GetPoint(idb,pointB)     
                    idbvero = point_locator.FindClosestPoint(pointB)
                    SolidpointEntityIdsArray.SetTuple1(idbvero,EdgeCellEntityId)
                    if j == (boundary.GetNumberOfPoints() - 1):
                        idb2 = int(boundariesExt.GetPointData().GetScalars().GetComponent(boundary.GetPointId(0),0))
                        meshToSurfaceExternal.Surface.GetPoint(idb2,pointB)
                        idb2vero = point_locator.FindClosestPoint(pointB)
                        EdgeCellArray2.InsertTuple3(index,idbvero+1, idb2vero+1, EdgeCellEntityId)
                        index = index + 1
                    else:
                        idb2 = int(boundariesExt.GetPointData().GetScalars().GetComponent(boundary.GetPointId(j+1),0))
                        meshToSurfaceExternal.Surface.GetPoint(idb2,pointB)
                        idb2vero = point_locator.FindClosestPoint(pointB)      
                        EdgeCellArray2.InsertTuple3(index,idbvero+1, idb2vero+1, EdgeCellEntityId)
                        index = index + 1
            #Write LifeV file with edges
            OutputFileName = self.LifeVSolidMeshName
            f=open(OutputFileName, 'w')
            line = "MeshVersionFormatted 1\n\n"
            line += "Dimension\n"
            line += "3\n\n"
            line += "Vertices\n"
            line += "%d\n" % Solid.GetNumberOfPoints()
            f.write(line)
            for i in range(Solid.GetNumberOfPoints()):
                point = Solid.GetPoint(i)
                pointValue = SolidpointEntityIdsArray.GetValue(i)
                line = "%f  %f  %f  %d\n" % (point[0]*self.ScaleFactor, point[1]*self.ScaleFactor, point[2]*self.ScaleFactor, pointValue)
                f.write(line)
            line = "\n"
            #tetras
            tetraCellIdArray = vtk.vtkIdTypeArray()
            Solid.GetIdsOfCellsOfType(tetraCellType,tetraCellIdArray)
            numberOfTetras = tetraCellIdArray.GetNumberOfTuples()
            line += "Tetrahedra\n"
            line += "%d\n" % numberOfTetras
            f.write(line)
            for i in range(numberOfTetras):
                tetraCellId = tetraCellIdArray.GetValue(i) 
                cellPointIds = Solid.GetCell(tetraCellId).GetPointIds()
                line = ''
                for j in range(cellPointIds.GetNumberOfIds()):
                    if j > 0:
                        line += '  '
                    line += "%d" % (cellPointIds.GetId(j)+1)
                cellEntityId = cellEntityIdsArray.GetValue(tetraCellId)   
                line += '  %d\n' % (cellEntityId+1)
                f.write(line)
            line = "\n"
            #tirangles
            triangleCellIdArray = vtk.vtkIdTypeArray()
            Solid.GetIdsOfCellsOfType(triangleCellType,triangleCellIdArray)
            numberOfTriangles = triangleCellIdArray.GetNumberOfTuples()
            line += "Triangles\n"
            line += "%d\n" % numberOfTriangles
            f.write(line)
            for i in range(numberOfTriangles):
                triangleCellId = triangleCellIdArray.GetValue(i)
                cellPointIds = Solid.GetCell(triangleCellId).GetPointIds()
                line = ''
                for j in range(cellPointIds.GetNumberOfIds()):
                    if j > 0:
                        line += '  '
                    line += "%d" % (cellPointIds.GetId(j)+1)
                cellEntityId = cellEntityIdsArray.GetValue(triangleCellId)
                line += '  %d\n' % cellEntityId
                f.write(line)
            #Edges  
            line += "Edges\n"
            line += "%d\n" % numberOfEdgesTot
            f.write(line)
            for i in range(EdgeCellArray1.GetNumberOfTuples()):
                line = "%d  %d  %d\n" % (EdgeCellArray1.GetComponent(i,0), EdgeCellArray1.GetComponent(i,1), EdgeCellArray1.GetComponent(i,2))
                f.write(line)
            for i in range(EdgeCellArray2.GetNumberOfTuples()):
                line = "%d  %d  %d\n" % (EdgeCellArray2.GetComponent(i,0), EdgeCellArray2.GetComponent(i,1), EdgeCellArray2.GetComponent(i,2))
                f.write(line) 
            f.close()

            print 'Solid LifeV mesh file written: ', self.LifeVSolidMeshName
            print ' '

            if self.VTKPrefixFilesName != None :
                nomeOutMeshVTU2 = "%s%s" % (self.VTKPrefixFilesName,"_SolidMesh.vtu")
                self.WriteUnstructuredGrid(Solid,nomeOutMeshVTU2)

            self.SolidMesh = Solid

if __name__=='__main__':
    main = pypes.pypeMain()
    main.Arguments = sys.argv
    main.Execute()      
