#Currently works only with version 5.10 of VTK. Several functions are used, which were removed in VTK 6.0

import vtk,ipdb
import numpy as np
import glob,os,sys,re
from vtk.util.numpy_support import vtk_to_numpy


def natural_sort(l):
    #needed to sort the filenames; taken from Stackoverflow
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    return sorted(l, key = alphanum_key)

class vtkFile:
    waterArray = {'Truchas':'VOF0002',
            'TruchasEnSight':'vof',
            'Thetis':'EAU',
            'Gerris':'T',
            'OpenFOAM':'alpha1'}
    velArray = {'Truchas':'Velocity',
            'TruchasEnSight':'velocity',
            'Thetis':'velocity',
            'Gerris':'Velocity',
            'OpenFOAM':'U'}
    pressureArray = {'Truchas':'P',
            'TruchasEnSight':'pressure',
            'Thetis':'pressure', #?
            'Gerris':'P',
            'OpenFOAM':'P'} #?
    
    def __init__(self,*args,**kwargs):
        self.solver=kwargs.get('solver','Truchas') #'Truchas' stands for default
        self.sim=kwargs.get('sim','2D') #'Truchas' stands for default
        if len(args)>0:
            path = args[0]
            if os.path.isfile(path) and os.path.exists(path):
                self.readFile(path)
    
    def readEnSight(self,ts=0,verbose=False):
        sys.exit('ReadEnSight function need improvement, access a connectionPort to a certain block of data')
        self.reader = vtk.vtkEnSight6Reader()
        self.reader.SetCaseFileName(self.filename)
        self.reader.ReadAllVariablesOn()
        self.reader.Update()
        timesteps = vtk_to_numpy(self.reader.GetTimeSets().GetItem(0))
        if ts >= len(timesteps):
            ts = len(timesteps) - 1 
        self.reader.SetTimeValue(timesteps[ts])
        #self.output = self.reader.GetOutput().GetBlock(0) #!!!
        #probably a filter choosing a single block of data is needed
    
    def readVTK(self,verbose=False):
        self.reader=vtk.vtkUnstructuredGridReader()
        try:
            self.reader.SetFileName(self.filename)
            self.reader.ReadAllScalarsOn()
            self.reader.ReadAllVectorsOn()
        except:
            print 'Could not find file: ',self.filename
            return
            
    def readFile(self,filename,verbose=False,ts = 0):  #ts required only for ensight files
        self.filename=os.path.abspath(filename)
        
        if os.path.splitext(self.filename)[1] == '.vtk':
            self.readVTK()
        elif os.path.splitext(self.filename)[1] == '.CASE':
            self.readEnSight(ts)
        self.reader.Update()
        self.output = self.reader.GetOutput()
        self.outputPort = self.reader.GetOutputPort()

        
        #a workaround: in binary vtk files (OpenFOAM, Truchas), the arrays are given, but with no information 
        #if they are scalars or vectors, therefore vectors are not rotated from 2D to 3D
        #it's probably a bug in current implementation of VTK
        self.pointData = self.output.GetPointData()
        self.cellData = self.output.GetCellData()
        if not self.pointData.GetVectors():
            self.pointData.SetVectors(self.pointData.GetArray(self.velArray[self.solver]))
            self.cellData.SetVectors(self.cellData.GetArray(self.velArray[self.solver]))
        
        #transformations:
        self.transform = vtk.vtkTransform() #transform is a matrix, and identity at the beginning, which is immediately multiplied by using Translate, Rotate etc.
        self.transform.PostMultiply()  #required to apply the transformation in the given order
        
        if self.sim=='2D':
            self.transform.RotateX(90)
            if self.solver in ['Gerris','Thetis']:
                self.transform.Scale(1.,0.,1.)
        if 'rotated' in self.filename:
            self.transform.Translate(0,0,-4.22)
            self.transform.RotateY(-10) #it is in XZ plane

        if self.solver == 'Gerris' and ('propagation' in self.filename or 'runup' in self.filename or 'stillwater' in self.filename):  
            self.transform.Translate(0,0,-1)#the water level is originally at z=1
        self.transformFilter=vtk.vtkTransformFilter()
        self.transformFilter.SetInputConnection(self.outputPort)
        self.transformFilter.SetTransform(self.transform)
        self.transformFilter.Update()
        self.output = self.transformFilter.GetOutput()
        self.outputPort = self.transformFilter.GetOutputPort()

        if self.output.GetCellData().GetNumberOfArrays() == 0:
            self.converter = vtk.vtkPointDataToCellData()
            self.converter.SetInputConnection(self.outputPort)
            self.converter.PassPointDataOn()
            self.converter.Update()
            self.output=self.converter.GetOutput() #later on output will always have at least cell data
            self.outputPort = self.converter.GetOutputPort()
        elif self.output.GetPointData().GetNumberOfArrays() == 0:
            self.converter = vtk.vtkCellDataToPointData()
            self.converter.SetInputConnection(self.outputPort)
            self.converter.PassCellDataOn()
            self.converter.Update()
            self.output=self.converter.GetOutput() #later on output will always have at least point data
            self.outputPort = self.converter.GetOutputPort()
        
        #self.output.Update()
        self.pointData = self.output.GetPointData()
        self.cellData  = self.output.GetCellData()
        
        
        self.getArrayNames(verbose)
        self.Ncells=self.output.GetNumberOfCells()
        
        #self.header = self.reader.GetHeader()
        self.time=-1

    def getArrayNames(self,verbose=False):
        self.arrNames=[self.cellData.GetArrayName(i) for i in range(self.cellData.GetNumberOfArrays())]
        if verbose:
            print 'Found',self.cellData.GetNumberOfArrays(),'arrays:',self.arrNames
        return self.arrNames

    def getCellArray(self,arrName):
        return vtk_to_numpy(self.cellData.GetArray(arrName))
        
    def getPointArray(self,arrName):
        return vtk_to_numpy(self.pointData.GetArray(arrName))
        
    def getVelocity(self):
        vel = self.getCellArray(self.velArray[self.solver])
        return vel
 
    def getWaterFraction(self):
        return self.getCellArray(self.waterArray[self.solver])
        
    def getPressure(self):
        return self.getCellArray(self.pressureArray[self.solver])
    
    def getCenters(self):
        vtkCenters=vtk.vtkCellCenters()
        vtkCenters.SetInputConnection(self.outputPort)
        vtkCenters.Update()
        centersOutput=vtkCenters.GetOutput()
        self.centers=np.array([centersOutput.GetPoint(i) for i in range(self.Ncells)])
        return self.centers
    
    def getVolumes(self):
        bounds=[self.output.GetCell(i).GetBounds() for i in range(self.Ncells)]
        self.bounds=np.array(bounds)
        self.dim=self.bounds[:,1::2]-self.bounds[:,0::2]  # ::2 means 'every second'
        if self.sim == '2D': 
            self.dim[:,1] = 1.
        self.vol= self.dim[:,0]*self.dim[:,1]*self.dim[:,2]

        return self.vol
        
    def getContour(self,value=0.5,arrName='',largestRegion=False):
        if arrName == '':
            arrName = self.waterArray[self.solver]
        contour=vtk.vtkContourFilter()
        contour.SetValue(0,value)
        self.pointData.SetActiveScalars(arrName)
        #contour.SetInputArrayToProcess(1, 0,0, vtkDataObject::FIELD_ASSOCIATION_POINTS , arrName);  #something like this may also work
        contour.SetInputConnection(self.outputPort)
        contour.Update()
        if largestRegion:
            conn=vtk.vtkConnectivityFilter()
            conn.SetInputConnection(contour.GetOutputPort())
            conn.SetExtractionModeToLargestRegion()
            conn.Update()
            connOutput = conn.GetOutput()
            IDs = []
            for iC in range(connOutput.GetNumberOfCells()):  #it has to be so complicated, no polylines are given
                cell = connOutput.GetCell(iC)
                for i in [0,1]:
                    ID = cell.GetPointId(i)
                    if not ID in IDs: IDs.append(ID)
            points=[connOutput.GetPoint(i) for i in IDs]
        else:
            cOutput = contour.GetOutput()
            points=[cOutput.GetPoint(i) for i in range(cOutput.GetNumberOfPoints())]
        
        if not points: sys.exit('vtkextract.py: No contours found - wrong initialization of water fraction?')
        
        points = np.array(points)
        return points
        
