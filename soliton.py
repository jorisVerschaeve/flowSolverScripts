import sys,os,glob,subprocess,re,shutil
import numpy as np
from scipy.interpolate import griddata
from ipdb import set_trace
from postprocess import readGphov as readGphov
import postprocess as pproc

#Grimshaw's formula for soliton speed
#c= sqrt( 1. + epsilon - 1./20.*(epsilon^2) - 3./70*(epsilon^3));

referenceDirectoryMain = os.path.join(os.getenv('HOME'),'work/soliton/fullPotentialSolution')
vofOfManyPolygons = os.path.join(os.getenv('HOME'),'work/soliton/initialize/bin/vofOfManyPolygons')

g = 9.80665
solitonVelFenton = {0.1:1.0485500*np.sqrt(g), 0.3:1.1376919*np.sqrt(g)}
solitonVelBIM =    {0.1:1.0485139*np.sqrt(g), 0.3:1.1375230*np.sqrt(g)}

def determineMesh():
    solver = pproc.determineSolver()
    if solver == 'Gerris':
        lines = open(glob.glob('*.gfs')[0]).readlines()
        sx = 1
        #FINISHED, BUT I COULD MOVE l = re.split(...) TO THE BEGINNING, TO MAKE IT CLEANER
        for l in lines:
            if 'GfsSimulation' in l:
                NBoxes = int(l.split()[0])
            if 'PhysicalParams' in l:
                H = float(re.split('[ =,{}]+',l.strip()) [2])
            if 'MetricStretch' in l and not '#' in l:
                sx = float(re.split('[ =,{}]+',l.strip()) [2])
            l = l.split()
            if len(l)==2 and l[0] in ['GfsRefine', 'Refine']  and l[1].isdigit():
                refine = int(l[1])
            if len(l)==3 and l[0] == 'Define' and l[1]== 'REFLEVEL':
                refine = int(l[2])
        Nx = NBoxes * 2**refine
        Ny = 2**refine
        L = H*sx*NBoxes
    elif solver == 'Truchas':
        lines = open(glob.glob('*.inp')[0]).readlines()
        for l in lines:
            l = filter(None, re.split('[ =,{}\n]+',l))
            if 'Ncell' in l:
                Nx = int(l[1])
                Ny = int(l[2])
            if 'Coord' in l:
                L = float(l[4]) - float(l[1])
                H = float(l[5]) - float(l[2])
    elif solver == 'Thetis':
        H = 0
        L = 0
        lines = open(glob.glob('*[0-9]*.don')[0]).readlines()
        for l in lines:
            l = re.split('[ ]+',l.strip())
            if len(l)==3 and l[0] == 'MAILLAGE':
                Nx = int(l[1])
                Ny = int(l[2])
            if 'DIM_MAX' in l:
                L = L + float(l[1].replace('D0',''))
                H = H + float(l[2].replace('D0',''))
            if 'DIM_MIN' in l:
                L = L - float(l[1].replace('D0',''))
                H = H - float(l[2].replace('D0',''))
    return L,H,Nx,Ny

def determineRotation():
    if 'rotated' in os.getcwd():
        if 'stillwater' in os.getcwd():
            angle = np.pi/18
            tx=0.5
            tz = -0.5
        else:
            a = pproc.determineAmplitude()
            angle = np.pi/18   #10 degrees rotation of the domain
            tx=0.0
            tz = 4.22
        return angle,tx,tz
    else:
        return 0,0,0

def rotateMesh(x,y,angle):
    X = x*np.cos(angle) - y*np.sin(angle)
    Y = x*np.sin(angle) + y*np.cos(angle)
    return X,Y

def meshCenters(Nx,Ny,x1,x2,y1,y2, tx=0,ty=0,angle=0):
    #tx,ty - translation of the coordinate system (before rotation, y - vertical direction)
    dx = (x2-x1)/Nx
    dy = (y2-y1)/Ny

    x,y=np.mgrid[x1+dx/2:x2-dx/2:Nx*1j,y1+dy/2:y2-dy/2:Ny*1j]
    x = x.transpose()
    y = y.transpose()
    
    #translating
    x = x - tx
    y = y - ty

    X,Y = rotateMesh(x,y,angle)
    return X,Y
    
def meshNodes(Nx,Ny,x1,x2,y1,y2, tx=0,ty=0,angle=0):
    #tx,ty - translation of the coordinate system (before rotation, y - vertical direction)
    dx = (x2-x1)/Nx
    dy = (y2-y1)/Ny

    x,y=np.mgrid[x1:x2:Nx*1j+1j,y1:y2:Ny*1j+1j]
    x = x.transpose()
    y = y.transpose()
    
    #translating
    x = x - tx
    y = y - ty

    X,Y = rotateMesh(x,y,angle)
    return X,Y
    

def initializeGerris(refDir,solitonPosition=15):
    L,H,Nx,Ny = determineMesh()
    dx = L/Nx
    dy = H/Ny
    angle,tx,ty = determineRotation()
    
    f = open('polygons.tmp','w')
    
    X,Y = meshNodes(Nx,Ny,0,L,-1,H-1,tx=tx,ty=ty,angle=angle)  #rotated mesh, but also translated to have interface at z=0
    eps = 0.000001  #required, so that interpolation points are close, but do not coincide with the actual nodes

    for i in range(X.shape[0]-1):
        for j in range(X.shape[1]-1):
            f.write( "{0:10.6f} {1:10.6f}".format(X[i  ,j  ],Y[i  ,j  ]))
            f.write(" {0:10.6f} {1:10.6f}".format(X[i+1,j  ],Y[i+1,j  ]))
            f.write(" {0:10.6f} {1:10.6f}".format(X[i+1,j+1],Y[i+1,j+1]))
            f.write(" {0:10.6f} {1:10.6f}".format(X[i,j+1  ],Y[i  ,j+1]))
            f.write("\n")
    f.close()
    
    subprocess.call(vofOfManyPolygons + ' interface.ini polygons.tmp t.cgd',shell=True)
    
    Xc,Yc = meshCenters(Nx,Ny,0-eps,L+eps,0-eps,H+eps) #actual mesh, interface is at z=1
    
    vof = open('t.cgd').read()
    f = open('t.cgd','w')
    f.write('2 y x\n')
    f.write('{0} {1}\n'.format(Ny,Nx))
    np.savetxt(f,Yc[:,0].reshape(1,-1),fmt='%12.6f') #in Gerris sw-level is at z=1; reshape to have it in 1 line
    np.savetxt(f,Xc[0,:].reshape(1,-1),fmt='%12.6f') #dx/2 needed to point at centers of cells
    f.write(vof)
    f.close()


    #VELOCITY
    Xc_rot, Yc_rot = meshCenters(Nx,Ny,0,L,-1,H-1,tx=tx,ty=ty,angle=angle)  #interface at z=0, for interpolation purposes
    cellCenters = np.array([Xc_rot.flatten(),Yc_rot.flatten()]).transpose()
    for component in ['u','v']:
        u,extent = readGphov(os.path.join(refDir,component))
        Ny_ref,Nx_ref = u.shape
        X_ref,Y_ref = meshNodes(Nx_ref-1,Ny_ref-1,*extent,tx=-solitonPosition)
        referencePoints = np.array([X_ref.flatten(),Y_ref.flatten()]).transpose()
        
        f = open(component+'.cgd','w')
        f.write('2 y x\n')
        f.write('{0} {1}\n'.format(Ny,Nx))
        interp = griddata(referencePoints, u.flatten(), cellCenters, fill_value=0)
        
        np.savetxt(f,Yc[:,0].reshape(1,-1),fmt='%12.6f') #in Gerris sw-level is at z=1; reshape to have it in 1 line
        np.savetxt(f,Xc[0,:].reshape(1,-1),fmt='%12.6f') #dx/2 needed to point at centers of cells
        np.savetxt(f,interp,fmt='%.4f') 
        f.close()
        

def initializeOpenFOAM(refDir,path='./',solitonPosition=15):
    points = open(os.path.join(path,'constant/polyMesh/points')).readlines()
    points = [p.replace('(','').replace(')','') for p in points if p[0] == '(']
    points = np.loadtxt(points)
    faces = open(os.path.join(path,'constant/polyMesh/faces')).readlines()
    #faces = [f[2:].replace(')','') for f in faces if f[:2] == '4(']
    faces = [f.replace(')','') for f in faces if f[:2] in ['3(','4(']]
    for i in range(len(faces)):
        if faces[i][:2] == '3(':
            faces[i] = faces[i][:-1] + ' ' + faces[i].split()[-1]+'\n'
        faces[i] = faces[i][2:]
    faces = np.loadtxt(faces,dtype=np.int)
    #reading 'owner' - assigning a cell number to each face
    owner = open(os.path.join(path,'constant/polyMesh/owner')).readlines()
    owner = [p for p in owner if p[0].isdigit()]
    owner = np.loadtxt(owner[1:],dtype=np.int)
    
    #3D->2D
    faces = points[faces.flatten()].reshape((-1,12)) #replacing indices with coordinates
    cellFilter=np.all(faces[:,2::3]==faces[:,2::3].max(),axis=1) #this is needed, because in OpenFOAM mesh is 3D, therefore 2D version must be obtained
    faces = faces.reshape((-1,3))[:,:2].reshape((-1,8))  #third dimension is not needed
    faces = faces[cellFilter]
    owner = owner[cellFilter]

    
    cells = np.zeros_like(faces)
    for iF,iC in enumerate(owner):
        cells[iC] = faces[iF]
        
    angle,tx,ty = determineRotation()
    x = cells[:,::2]
    y = cells[:,1::2]
    x = x - tx
    y = y - ty
    cells[:,::2],cells[:,1::2] = rotateMesh(x,y,angle)
    
    np.savetxt('polygons.tmp',cells,fmt='%12.6f')
    subprocess.call(vofOfManyPolygons + ' interface.ini polygons.tmp vof.ini',shell=True)

    #formatting alpha1 OpenFOAM initial field file
    org = open(os.path.join(path,'0/alpha1')).readlines()
    f = open(os.path.join(path,'0/alpha1'),'w')
    skip=0
    for i,l in enumerate(org):
        if not skip:
            if 'internalField' in l:
                if l.split()[1] == 'nonuniform':
                    f.write(l)
                    skip = int(org[i+1]) + 4   #counts for old values, opening and closing brackets and the number of cells and a semicolon
                else:
                    #replacing the 'internalField   uniform 0;' line
                    f.write('internalField   nonuniform List<scalar>\n')

                f.write(str(len(cells)) + '\n(\n')
                f.write(open('vof.ini').read())
                f.write(')\n;\n')
            else:
                f.write(l)
        else:
            skip = skip-1
    f.close()

    # VELOCITY AND PRESSURE
    interp = {}
    cellCenters = np.array([np.average(cells[:,::2],axis=1),np.average(cells[:,1::2],axis=1)]).transpose()
    for component in ['u','v','p']:
        u,extent = readGphov(os.path.join(refDir,component),rescale=False)
        Ny_ref,Nx_ref = u.shape
        X_ref,Y_ref = meshNodes(Nx_ref-1,Ny_ref-1,*extent,tx=-solitonPosition)
        referencePoints = np.array([X_ref.flatten(),Y_ref.flatten()]).transpose()
        interp[component] = griddata(referencePoints, u.flatten(), cellCenters, fill_value=0)
        
    interp['p'] = (interp['p'] - cellCenters[:,1]) * g*1000
    interp['u'] =  interp['u'] * np.sqrt(g)
    interp['v'] =  interp['v'] * np.sqrt(g)
    org = open(os.path.join(path,'0/U')).readlines()
    f = open(os.path.join(path,'0/U'),'w')
    skip=0
    for i,l in enumerate(org):
        if not skip:
            if 'internalField' in l:
                if l.split()[1] == 'nonuniform':
                    f.write(l)
                    skip = int(org[i+1]) + 4   #counts for old values, opening and closing brackets and the number of cells and a semicolon
                else:
                    #replacing the 'internalField   uniform 0;' line
                    f.write('internalField   nonuniform List<vector>\n')

                f.write(str(len(cells)) + '\n(\n')
                np.savetxt(f,np.array([interp['u'],interp['v'],np.zeros(len(cells))]).transpose(),fmt='(%10.6f %10.6f %10.6f)')
                f.write(')\n;\n')
            else:
                f.write(l)
        else:
            skip = skip-1
    f.close()

    org = open(os.path.join(path,'0/p_rgh')).readlines()
    f = open(os.path.join(path,'0/p_rgh'),'w')
    skip=0
    for i,l in enumerate(org):
        if not skip:
            if 'internalField' in l:
                if l.split()[1] == 'nonuniform':
                    f.write(l)
                    skip = int(org[i+1]) + 4   #counts for old values, opening and closing brackets and the number of cells and a semicolon
                else:
                    #replacing the 'internalField   uniform 0;' line
                    f.write('internalField   nonuniform List<scalar>\n')

                f.write(str(len(cells)) + '\n(\n')
                np.savetxt(f,interp['p'],fmt='%10.6f')
                f.write(')\n;\n')
            else:
                f.write(l)
        else:
            skip = skip-1
    f.close()

def initializeTruchas(refDir, path='./', solitonPosition=15):
    L,H,Nx,Ny = determineMesh()
    angle,tx,ty = determineRotation()
    X,Y = meshNodes(Nx,Ny,0,L,-1,H-1,tx=tx,ty=ty,angle=angle)
    f = open('polygons.tmp','w')
    for i in range(X.shape[0]-1):
        for j in range(X.shape[1]-1):
            f.write( "{0:10.6f} {1:10.6f}".format(X[i  ,j  ],Y[i  ,j  ]))
            f.write(" {0:10.6f} {1:10.6f}".format(X[i+1,j  ],Y[i+1,j  ]))
            f.write(" {0:10.6f} {1:10.6f}".format(X[i+1,j+1],Y[i+1,j+1]))
            f.write(" {0:10.6f} {1:10.6f}".format(X[i,j+1  ],Y[i  ,j+1]))
            f.write("\n")
    f.close()

    subprocess.call(vofOfManyPolygons + ' interface.ini polygons.tmp vof.ini',shell=True)

    # VELOCITY
    Xc,Yc = meshCenters(Nx,Ny,0,L,-1,H-1,tx=tx,ty=ty,angle=angle)
    cellCenters = np.array([Xc.flatten(),Yc.flatten()]).transpose()
    for component in ['u','v']:
        u,extent = readGphov(os.path.join(refDir,component))
        Ny_ref,Nx_ref = u.shape
        X_ref,Y_ref = meshNodes(Nx_ref-1,Ny_ref-1,*extent,tx=-solitonPosition)
        referencePoints = np.array([X_ref.flatten(),Y_ref.flatten()]).transpose()
        interp = griddata(referencePoints, u.flatten(), cellCenters, fill_value=0)
        np.savetxt(component+'.ini',interp,fmt='%10.6f')

def initializeThetis(refDir, path='./', solitonPosition=15):
    L,H,Nx,Ny = determineMesh()
    angle,tx,ty = determineRotation()
    dx = L/Nx
    dy = H/Ny
    #Thetis needs values given for nodes,not cell centers; 
    #I have to contruct the volumes around 'nodes', to compute volume fractions for them
   
    os.mkdir('objets')
    interface = np.loadtxt('interface.ini')
    interface = np.vstack((interface[0],interface)) 
    interface = np.vstack((interface,interface[-1]))
    interface[0,1] = -20  #to create a closed polygon
    interface[-1,1] = -20
    interface[:,0],interface[:,1] = rotateMesh(interface[:,0],interface[:,1],-angle)
    interface[:,1] = interface[:,1]+ty
    N_intf = len(interface)
    f = open('objets/interface.mxa','w')
    f.write('{0} {0}\n'.format(N_intf))
    np.savetxt(f,interface,fmt='%10.4f')
    l = np.arange(1,N_intf+1)
    lista = np.vstack((l[:-1],l[1:])).transpose()
    np.savetxt(f,lista,fmt='%d')
    f.write('{0} {1}\n'.format(N_intf,1))
    f.close()
    
    Xc,Yc = meshNodes(Nx,Ny,0,L,-1,H-1,tx=tx,ty=ty,angle=angle)
    nodes = np.array([Xc.flatten(),Yc.flatten()]).transpose()
    for component in ['u','v']:
        u,extent = readGphov(os.path.join(refDir,component))
        Ny_ref,Nx_ref = u.shape
        X_ref,Y_ref = meshNodes(Nx_ref-1,Ny_ref-1,*extent,tx=-solitonPosition)
        referencePoints = np.array([X_ref.flatten(),Y_ref.flatten()]).transpose()
        interp = griddata(referencePoints, u.flatten(), nodes, fill_value=0)
        np.savetxt(component+'.ini',interp,fmt='%10.6f')
    
def plotSolitonIntersection():
    """Script plots several cells given in cells.txt file and a soliton surface crossing them.
    Used only for testing purposes, not needed any more."""
    
    cells = np.loadtxt('cells.txt')   #file with several cell [x0 y0 x1 y1 x2 y2 x3 y3], each cell in seperate row (at least 2 needed)
    cells = np.hstack([cells,cells[:,:2]])
    for cell in cells:
        plt.plot(cell[::2],cell[1::2])

    np.savetxt('X.txt',sorted(cells.flatten()[::2]))
    subprocess.call('solitonEta 0.3 15 0 X.txt',shell=True)
    soliton = np.loadtxt('eta.txt')
    plt.plot(soliton[:,0],soliton[:,1])

    #ipdb.set_trace()
    plt.savefig('boxes.png')
    plt.show()

    
def initialize(a=None, solitonPosition=15):
    print 'Initialization of a solitary wave from numerical data'
    sys.stdout.flush()
    solver = pproc.determineSolver()
    a = a or pproc.determineAmplitude()
    refDir = os.path.join(referenceDirectoryMain,str(a))
    if not os.path.exists(refDir):
        print 'No data for given amplitude found. Amplitude provided:',a,', available data:'
        for d in sorted(glob.glob(os.path.join(refDir,'../0.*'))):
            print d
        sys.exit()
    surf = np.loadtxt(os.path.join(refDir,'surf'))
    surf[:,0] = surf[:,0] + solitonPosition
    surf = np.vstack((np.array([-200,0.]),surf,np.array([200,0]))) #if the left side is inside the domain, vofOfManyPolygons has problems
    np.savetxt('interface.ini',surf,fmt='%10.6f')
        
    if solver == 'Gerris':
        initializeGerris(refDir)
    elif solver == 'OpenFOAM':
        initializeOpenFOAM(refDir)
    elif solver == 'Truchas':
        initializeTruchas(refDir)
    elif solver == 'Thetis':
        initializeThetis(refDir)

    
