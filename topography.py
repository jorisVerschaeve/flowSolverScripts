import sys,os,glob,subprocess,string, matplotlib,ipdb
import numpy as np
import matplotlib.pyplot as plt

#Converts the xyz topography file to formats readable by ParaView. 
#Currently works only for .csv format


def xyz2csv(path):
    """To open the csv file in ParaView, load it, Filters->Alphabetical->Table to Points
    ->Choose colums for x, y, z coordinates
    ->Split the view, select View 3D, make the filter visible (required due to a bug)"""
    xyz = np.loadtxt(path)
    print 'Succesfully loaded',path
    #Np = 1000000
    Np = len(xyz)
    print 'Writing to .csv'
    xyz = np.hstack([xyz,xyz[:,2].reshape((-1,1))])
    f = open(os.path.splitext(path)[0]+'.csv','w')
    f.write("x coord,y coord,z coord,height\n")
    #np.savetxt(f,xyz[:Np],fmt='%.2f,%.2f,%.2f,%.2f')
    np.savetxt(f,xyz[:Np],fmt='%.1f,%.1f,%d,%d')
    f.close()
    
    #cuts the domain to certain size
    xyz = xyz[xyz[:,0]>386000] 
    xyz = xyz[xyz[:,0]<420000]
    xyz = xyz[xyz[:,1]<6913000]
    f = open(os.path.splitext(path)[0]+'_mini.csv','w')
    f.write("x coord,y coord,z coord,height\n")
    #np.savetxt(f,xyz[:Np],fmt='%.2f,%.2f,%.2f,%.2f')
    np.savetxt(f,xyz[:Np],fmt='%.1f,%.1f,%d,%d')
    f.close()
    
    
    
def xyz2ply(xyzPath):  #not really working, would require faces to visualize in ParaView
    xyz = np.loadtxt(xyzPath)
    print 'Succesfully loaded',xyzPath
    ply = open(os.path.splitext(xyzPath)[0]+'.ply','w')
    #Np = len(xyz)
    Np = 100
    Nf = 0
    ply.write("""ply
format ascii 1.0
element vertex {0}
property float x
property float y
property float z
element face {1}
property list uchar int vertex_index
end_header
""".format(Np,Nf))
    np.savetxt(ply,xyz[:Np],fmt='%.4f %.4f %.4f')
    ply.close()
