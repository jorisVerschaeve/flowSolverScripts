import sys,os,glob,subprocess,re,shutil,ipdb
import numpy as np
from scipy.interpolate import griddata


refDir = os.path.join(os.getenv("HOME"),'work/soliton/slide_0.5_5_6_5_4/nf16/')

foamHeader="""/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  2.1.1                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
"""
header = foamHeader + """    class       vectorAverageField;
    object      values;
}

// Average
(0 0 0)

// Data on points
"""

headerPoints = foamHeader + """    class       vectorField;
    object      points;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
"""


def gaussian0(caseDir='./',xc=3):
    #creates gaussian data for OpenFOAM
    caseDir = os.path.abspath(caseDir)
    caseDir = os.path.join(caseDir,'constant/polyMesh')
    mesh = open(os.path.join(caseDir,'blockMeshDict')).readlines()
    for i,l in enumerate(mesh):
        l = l.split()
        if l:
            if l[0] == 'x1': x1 = float(l[1][:-1])
            if l[0] == 'x2': x2 = float(l[1][:-1])

    N = 100
    A = np.zeros((N,3))
    A[:,0] = np.linspace(x1,x2,N)
    A = A[1:-1,:]
    c = 0.2
    A[:,1] = 0.5*np.exp(-(A[:,0]-xc)**2/c**2/2)

    def save(name,A):
        np.savetxt(name,A,fmt='(%.6f %.6f %.1f)')

    save(os.path.join(caseDir,'gaussian1'),A)
    A[:,2] = 1
    save(os.path.join(caseDir,'gaussian2'),A)



def gaussianTimeSeries(caseDir = './',xc=3):
    caseDir = os.path.abspath(caseDir)
    
    gaussian0(caseDir,xc)
    
    T = np.linspace(0,1,21)
    dir = os.path.join(caseDir,'constant/boundaryData/lowerWall')
    subprocess.call('mkdir -p '+dir,shell=True)
    N = 200
    x1 = -1.5
    x2 = 4.5
    points = np.zeros((N,3))
    points[:,0] = np.linspace(x1,x2,N)
    points = np.vstack((points,points))
    points[N:,2] = 1
    
    f = open(os.path.join(dir,'points'),'w')
    f.write(headerPoints)
    f.write('(\n')
    np.savetxt(f,points,fmt='(%.6f %.6f %.1f)')
    f.write(')\n')
    f.close()
    
    #Gaussian parameters:
    c = 0.2
    
    #gauss0 is the initial shape of the boundary
    gauss0 = 0.5*np.exp(-(points[:,0]-xc)**2/c**2/2)
    
    for t in T:
        tDir = os.path.join(dir,"{:g}".format(t))
        subprocess.call('mkdir -p '+tDir,shell=True)
        f = open(os.path.join(tDir,'pointDisplacement'),'w')
        f.write(header)
        
        #should be interpolated eventually
        xct = xc - 0.5*t
        
        if t>=0:
            gauss = 0.5*np.exp(-(points[:,0]-xct)**2/c**2/2)   # - gauss0
        else:
            gauss = np.zeros(2*N)
        #zero on the edges of the domain:
        gauss[0] = 0
        gauss[N-1] = 0
        gauss[N] = 0
        gauss[-1] = 0
        f.write('(\n')
        np.savetxt(f,gauss,fmt='(0 %.6f 0)')
        f.write(')\n')
        f.close()
        
def timeSeriesFromReference(caseDir = './'):
    caseDir = os.path.abspath(caseDir)
    
    loc = np.loadtxt(os.path.join(refDir,'sb.dat'))
    negativeN = np.sum(loc[:,0]<=0)
    loc = loc[negativeN-1:,:]  #allow a zero or a negative time
    if loc[0,0] != 0:
        loc[0,0] = 0 #to guarantee a start from t=0
    
    #rescalling, as the location of the body is given in the nondimensional time:
    loc[:,0] = loc[:,0]*np.sqrt(1/9.81)
    disp = np.loadtxt(os.path.join(refDir,'sr.dat'))
    
    dir = os.path.join(caseDir,'constant/boundaryData/lowerWall')
    subprocess.call('mkdir -p '+dir,shell=True)
    N = 200
    mesh = open(os.path.join(caseDir,'constant/polyMesh/blockMeshDict')).readlines()
    for i,l in enumerate(mesh):
        l = l.split()
        if l:
            if l[0] == 'x1': x1 = float(l[1][:-1])
            if l[0] == 'x2': x2 = float(l[1][:-1])
    points = np.zeros((N,3))
    points[:,0] = np.linspace(x1,x2,N)
    points = np.vstack((points,points))
    points[N:,2] = 1
    
    f = open(os.path.join(dir,'points'),'w')
    f.write(headerPoints)
    f.write('(\n')
    np.savetxt(f,points,fmt='(%.6f %.6f %.1f)')
    f.write(')\n')
    f.close()
    
    
    for t,dx in loc:
        tDir = os.path.join(dir,"{:g}".format(t))
        subprocess.call('mkdir -p '+tDir,shell=True)
        f = open(os.path.join(tDir,'pointDisplacement'),'w')
        f.write(header)
        
        dispInterp = griddata(disp[:,0] + dx, disp[:,1], points[:,0], fill_value=0)
        
        #zeros on the edges of the domain:
        dispInterp[0] = 0
        dispInterp[N-1] = 0
        dispInterp[N] = 0
        dispInterp[-1] = 0
        f.write('(\n')
        np.savetxt(f,dispInterp,fmt='(0 %.6f 0)')
        f.write(')\n')
        f.close()
    
    
def slideHorizontal(caseDir = './'):
    from scipy.interpolate import interp1d
    caseDir = os.path.abspath(caseDir)
    angle=np.pi/4
    Nt=10
    loc = np.zeros((Nt,2))

    if loc[0,0] != 0:
        loc[0,0] = 0 #to guarantee a start from t=0
    
    disp = np.loadtxt(os.path.join(refDir,'sr.dat'))
    
    dir = os.path.join(caseDir,'constant/boundaryData/lowerWall')
    subprocess.call('mkdir -p '+dir,shell=True)
    mesh = open(os.path.join(caseDir,'constant/polyMesh/blockMeshDict')).readlines()
    xc = 0.5
    for i,l in enumerate(mesh):
        l = l.split()
        if l:
            if l[0] == 'x1': x1 = float(l[1][:-1])
            if l[0] == 'y1': y1 = float(l[1][:-1])
            if l[0] == 'x2': x2 = float(l[1][:-1])
            if l[0] == 'xc': xc = float(l[1][:-1])
    N=500
    points = np.zeros((N,3))
    points[:N/2,0] = np.linspace(0,1.5,N/2+1)[:-1]
    points[:N/2,1] = np.linspace(1.5,0,N/2+1)[:-1]
    points[N/2:,0] = np.linspace(1.5,3.5,N/2)
    points = np.vstack((points,points))
    points[N:,2] = 1
    Nsh = 5
    loc = np.zeros((Nsh,2))
    loc[:,0] = [0.0,0.1,0.3,0.4,1.5]
    loc[:,1] = [0.0,1.1,1.3,1.4,1.5]
    disp = np.zeros((Nsh,3))
    disp[:] = [[0.0,0.0,0.0],[0.0,0.0,0.0],[0.00,0.07,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]]
    t0=0
    t1=3
    T = np.linspace(t0,t1,151)
    #Interpolate the location of the bump, between two points of time
    interpLoc = interp1d([t0,t1],loc)
    f = open(os.path.join(dir,'points'),'w')
    f.write(headerPoints)
    f.write('(\n')
    np.savetxt(f,points,fmt='(%.6f %.6f %.6f)')
    f.write(')\n')
    f.close()
    distance = np.sqrt(points[:,0]**2+(points[:,1]-y1)**2)
    for t in T:
        tDir = os.path.join(dir,"{:g}".format(t))
        subprocess.call('mkdir -p '+tDir,shell=True)
        f = open(os.path.join(tDir,'pointDisplacement'),'w')
        f.write(header)
        
        dispInterp = griddata(interpLoc(t),disp,distance, fill_value=0)
        #dispInterp = np.zeros((400,1))
        
        #zeros on the edges of the domain:
        f.write('(\n')
        np.savetxt(f,dispInterp,fmt='(%.6f %.6f %.6f)')
        f.write(')\n')
        f.close()
    
    
def slideRotated(caseDir = './'):
    from scipy.interpolate import interp1d
    print "Generating point displacement time series"
    caseDir = os.path.abspath(caseDir)
    angle=np.pi/4
    Nt=10
    loc = np.zeros((Nt,2))

    dir = os.path.join(caseDir,'constant/boundaryData/lowerWall')
    subprocess.call('mkdir -p '+dir,shell=True)
    mesh = open(os.path.join(caseDir,'constant/polyMesh/blockMeshDict')).readlines()
    xc = 0.5
    for i,l in enumerate(mesh):
        l = l.split()
        if l:
            if l[0] == 'x1': x1 = float(l[1][:-1])
            if l[0] == 'y1': y1 = float(l[1][:-1])
            if l[0] == 'x2': x2 = float(l[1][:-1])
            if l[0] == 'xc': xc = float(l[1][:-1])
    N=500
    points = np.zeros((N,3))
    points[:,0] = np.linspace(0,1.9,N) #do not have to cover whole patch
    points = np.vstack((points,points))
    points[N:,2] = 1
    Nsh = 6
    loc = np.zeros((Nsh,2))
    loc[:,0] = [0.0,0.1,0.2,0.4,0.5,100]
    loc[:,1] = loc[:,0]
    loc[1:-1,1] += 1.4
    disp = np.zeros((Nsh,3))
    A = 0.15
    disp[:,1] = [0.0,0.0,A,A,0.0,0.0]
    t0=0
    t1=1.5
    t_sim = 5
    T = np.linspace(t0,t1,51)
    #Interpolate the location of the bump, between two points of time
    interpLoc = interp1d([t0,t1],loc)
    f = open(os.path.join(dir,'points'),'w')
    f.write(headerPoints)
    f.write('(\n')
    np.savetxt(f,points,fmt='(%.6f %.6f %.6f)')
    f.write(')\n')
    f.close()
    for t in T:
        tDir = os.path.join(dir,"{:g}".format(t))
        subprocess.call('mkdir -p '+tDir,shell=True)
        f = open(os.path.join(tDir,'pointDisplacement'),'w')
        f.write(header)
        
        dispInterp = griddata(interpLoc(t),disp,points[:,0], fill_value=0)
        #dispInterp = np.zeros_like(dispInterp)
        
        #zeros on the edges of the domain:
        f.write('(\n')
        np.savetxt(f,dispInterp,fmt='(%.6f %.6f %.6f)')
        f.write(')\n')
        f.close()
    subprocess.call('cp -r {} {}'.format(tDir,os.path.join(dir,"{:g}".format(t_sim))),shell=True)
   



vel_exp={1:2.45, 2:3.38,3:3.56}
h_exp = {1: 0.16, 2: 0.16, 3: 0.12}

def slideExpRotated(caseDir = './',exp=1):
    from scipy.interpolate import interp1d
    print "Generating point displacement time series"
    caseDir = os.path.abspath(caseDir)
    Nt=10
    loc = np.zeros((Nt,2))

    dir = os.path.join(caseDir,'constant/boundaryData/lowerWall')
    subprocess.call('mkdir -p '+dir,shell=True)
    mesh = open(os.path.join(caseDir,'constant/polyMesh/blockMeshDict')).readlines()
    xc = 0.5
    for i,l in enumerate(mesh):
        l = l.split()
        if l:
            if l[0] == 'x1': x1 = float(l[1][:-1])
            if l[0] == 'y1': y1 = float(l[1][:-1])
            if l[0] == 'x2': x2 = float(l[1][:-1])
            if l[0] == 'xc': xc = float(l[1][:-1])
    N=2000
    points = np.zeros((N,3))
    points[:,0] = np.linspace(0,1.4,N) #do not have to cover whole patch
    points = np.vstack((points,points))
    points[N:,2] = 1
    Nsh = 6
    loc = np.zeros((Nsh,2))
    h = h_exp[exp]
    vel = vel_exp[exp]
    x_back = -10
    x_front = 0.2
    dist = 1.2 #distance measured along the slope
    loc[:,0] = [-100,x_back,x_back+h,x_front-h,x_front,100] #x_front - h gives 45 degrees
    loc[:,1] = loc[:,0]
    loc[1:-1,1] += dist
    disp = np.zeros((Nsh,3))
    disp[:,1] = [0.0,0.0,h,h,0.0,0.0]
    alpha = 35*np.pi/180
    t1=dist/vel   
    T = np.linspace(0,t1,501)
    #Interpolate the location of the bump, between two points of time
    interpLoc = interp1d([0,t1],loc)
    f = open(os.path.join(dir,'points'),'w')
    f.write(headerPoints)
    f.write('(\n')
    np.savetxt(f,points,fmt='(%.6f %.6f %.6f)')
    f.write(')\n')
    f.close()
    for t in T:
        tDir = os.path.join(dir,"{:g}".format(t))
        subprocess.call('mkdir -p '+tDir,shell=True)
        f = open(os.path.join(tDir,'pointDisplacement'),'w')
        f.write(header)
        
        dispInterp = griddata(interpLoc(t),disp,points[:,0], fill_value=0)
        #dispInterp = np.zeros_like(dispInterp)
        
        #zeros on the edges of the domain:
        f.write('(\n')
        np.savetxt(f,dispInterp,fmt='(%.6f %.6f %.6f)')
        f.write(')\n')
        f.close()
    
def slideExpHorizontal(caseDir = './',exp=1):
    from scipy.interpolate import interp1d
    caseDir = os.path.abspath(caseDir)
    dir = os.path.join(caseDir,'constant/boundaryData/lowerWall')
    subprocess.call('mkdir -p '+dir,shell=True)
    mesh = open(os.path.join(caseDir,'constant/polyMesh/blockMeshDict')).readlines()
    xc = 0.5
    for i,l in enumerate(mesh):
        l = l.split()
        if l:
            if l[0] == 'x1': x1 = float(l[1][:-1])
            if l[0] == 'y1': y1 = float(l[1][:-1])
            if l[0] == 'x2': x2 = float(l[1][:-1])
            if l[0] == 'xc': xc = float(l[1][:-1])
    N=500
    points = np.zeros((N,3))
    points[:N/2,0] = np.linspace(0,xc,N/2+1)[:-1]
    points[:N/2,1] = np.linspace(y1,0,N/2+1)[:-1]
    points[N/2:,0] = np.linspace(xc,3.5,N/2)
    points = np.vstack((points,points))
    points[N:,2] = 1
    Nsh = 6
    x_front = 1.0
    x_back = 0.1
    h = 0.12
    loc = np.zeros((Nsh,2))
    alpha=35*np.pi/180
    loc[:,0] = [0.0,x_back,x_back+h*np.cos(alpha) + h*np.sin(alpha),x_front - h*np.cos(alpha) + h*np.sin(alpha),x_front,1.5]
    dist = 0.5
    vel = vel_exp[exp]
    loc[:,1] = loc[:,0]
    loc[1:-1,1] += dist * np.cos(alpha)

    disp = np.zeros((Nsh,3))
    #disp[:,0] = [0.0,0.0,sin*h,sin*h,0.0,0.0]
    disp[:,1] = [0.1,0.1,h/np.cos(alpha),h/np.cos(alpha),0.0,0.0]
    t0=0
    t1=dist/vel
    T = np.linspace(t0,t1,51)
    #Interpolate the location of the bump, between two points of time
    interpLoc = interp1d([t0,t1],loc)
    f = open(os.path.join(dir,'points'),'w')
    f.write(headerPoints)
    f.write('(\n')
    np.savetxt(f,points,fmt='(%.6f %.6f %.6f)')
    f.write(')\n')
    f.close()
    distance = np.sqrt(points[:,0]**2+(points[:,1]-y1)**2)
    for t in T:
        tDir = os.path.join(dir,"{:g}".format(t))
        subprocess.call('mkdir -p '+tDir,shell=True)
        f = open(os.path.join(tDir,'pointDisplacement'),'w')
        f.write(header)
        
        dispInterp = griddata(interpLoc(t),disp,distance, fill_value=0)
        #dispInterp = np.zeros((400,1))
        
        #zeros on the edges of the domain:
        f.write('(\n')
        np.savetxt(f,dispInterp,fmt='(%.6f %.6f %.6f)')
        f.write(')\n')
        f.close()
    
    
