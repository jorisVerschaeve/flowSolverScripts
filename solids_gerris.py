import numpy as np
import subprocess,sys,os
from scipy.interpolate import griddata

from solids_foam import refDir

def cylinder(fn='cylinder'):  #the body given as z = f(x,y)
    def cylinderFunction(x,y,R):
        res = np.zeros(len(x))
        fi = np.arccos(x[np.abs(x)<=R]/R)
        res[np.abs(x)>R] = -np.abs(x/R)
        res[np.abs(x)<=R] = R*np.sin(fi)
        return res
    N=100
    M = 10
    R = 1
    A = np.zeros((N*M,3))
    Rm=1.1*R
    X,Y = np.mgrid[-Rm:Rm:N*1j,-0.4:0.4:M*1j]
    A[:,0] = X.flatten()
    A[:,1] = Y.flatten()
    A[:,2] = cylinderFunction(A[:,0],A[:,1],R)

    f = open(fn+'.xyz','w')
    f.write('{} 0 0\n'.format(N*M))
    np.savetxt(f,A,fmt='%.6f %.6f %.6f')
    f.close()
    



def triangulate(fn='cylinder'):
    subprocess.call('cat {0}.xyz | delaunay > {0}.gts'.format(fn),shell=True)
    subprocess.call('gtscheck -v < {}.gts'.format(fn),shell=True)

    # extract a contour level at z = 0
    subprocess.call('gts2oogl -G -n -I 0 < {}.gts > iso'.format(fn),shell=True)
    # add the points of this contour to the original dataset and retriangulate
    iso = open('iso').readlines()
    iso = ''.join([i for i in iso if len(i)>0] )
    cylinder = open(fn+'.xyz').readlines()
    N = int(cylinder[0].split()[0])
    cylinder[0]  = cylinder[0].replace(str(N),str(N+len(np.loadtxt('iso')) ) )
    cylinder = ''.join([i for i in cylinder if len(i)>0] )
    f = open(fn+'.xyz','w')
    f.write(cylinder)
    f.write(iso)
    f.close()
    subprocess.call("cat {0}.xyz | delaunay -d > {0}.gts".format(fn),shell=True)
    subprocess.call('gtscheck -v < {}.gts'.format(fn),shell=True)


def timeSeriesFromReference(caseDir = './',fn='solid'): #the body given as z = f(x,y)
    caseDir = os.path.abspath(caseDir)
    
    #MOVEMENT OF THE BODY
    #the location of the file will be given as a Cartesian Grid Data (CGD) file (look for documentation of GfsFunction)
    loc = np.loadtxt(os.path.join(refDir,'sb.dat'))
    #rescalling, as the location of the body is given in the nondimensional time:
    loc[:,0] = loc[:,0]*np.sqrt(1/9.81)
    dx0 = griddata(loc[:,0],loc[:,1],0)  #initial translation (to be added to the shape of the body)

    
    #Gerris requires the velocity of the moving surface, not the location
    dt = np.diff(loc[:,0])
    vel = np.diff(loc[:,1])/dt
    tVel = loc[:-1,0] + dt/2
    
    #Gerris requires the data be given on a regular mesh (with constant dt), so we interpolate
    M = 101
    T1 = 10.0
    dt = T1/M
    T = np.linspace(0,T1,M)
    
    velInterp = griddata(tVel, vel, T, fill_value=0)
    
    f = open('velocity.cgd','w')
    f.write('1 t\n{}\n'.format(M))
    np.savetxt(f,T.reshape((1,-1)),fmt=' %.2f')
    np.savetxt(f,velInterp,fmt=' %.4f')
    f.close()
    
    
    #SHAPE OF THE BODY
    points = np.loadtxt(os.path.join(refDir,'sr.dat'))
    points[:,0] += dx0
    N = len(points)
    points = np.transpose(np.vstack((points[:,0],np.zeros(N),points[:,1])))
    points = np.vstack((np.array([np.min(points[:,0])-0.1,0,-0.1]),points,np.array([np.max(points[:,0])+0.1,0,-0.1])))
    N = N+2
    points = np.vstack((points,points))
    points[:N,1] = -1
    points[N:,1] = 1
    
    f = open(os.path.join(caseDir,fn+'.xyz'),'w')
    f.write('{} 0 0\n'.format(N*2))
    np.savetxt(f,points,fmt='%.6f %.6f %.6f')
    f.close()
 
vel_exp={1:2.45, 2:3.38,3:3.56}
h_exp = {1: 0.16, 2: 0.16, 3: 0.12}

   
    
def slideExp(caseDir = './',exp=1): 
    caseDir = os.path.abspath(caseDir)
    
    h = h_exp[exp]
    vel = vel_exp[exp]
    x_back = -10
    x_front = 0.2
    dist = 1.2 #distance measured along the slope
    N = 6
    points = np.zeros((N,3))
    points[:,0] = [x_back-0.1,x_back,x_back+h,x_front-h,x_front,x_front+0.1] #x_front - h gives 45 degrees
    points[:,2] = [-0.1,0.0,h,h,0.0,-0.1]
    points = np.vstack((points,points))
    points[:N,1] = -1
    points[N:,1] = 1

    t1=dist/vel
    T = np.linspace(0,t1,51)

    f = open(os.path.join(caseDir,'body.xyz'),'w')
    f.write('{} 0 0\n'.format(N*2))
    np.savetxt(f,points,fmt='%.6f %.6f %.6f')
    f.close()
    
    
    
