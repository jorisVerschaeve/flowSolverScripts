#!/usr/bin/env python

import sys,os, glob, subprocess, string, matplotlib, ipdb
#matplotlib.use('pdf')
import postprocess
from phd import vtkextract
from phd.soliton import solitonVelBIM, solitonVelFenton, g
from ipdb import set_trace
import numpy as np
from scipy.interpolate import griddata
#matplotlib.use('Agg')

import matplotlib.pyplot as plt
from itertools import *
import vtk
matplotlib.rc('text', usetex=True)

solitonVel = solitonVelBIM


def sample_vtk(file, p1=None, p2=None, x_c=0, N_points=101, solver=None, verbose=True):
    """
    Samples velocity field in a given VTK file, along a line defined with p1 and p2, or x_c.
    :param file: Path to the VTK file
    :type file: str
    :param p1: Coordinates of the starting point of the sample line. Can be either 2D or 3D [x, y] or [x, y, z]
    :type p1: list
    :param p2: Coordinates of the ending point of the sample line. Can be either 2D or 3D [x, y] or [x, y, z]
    :type p2: list
    :param x_c: If p1 and p2 are not provided, this parameter serves the x coordinate of the vertical line.
    :type x_c: float
    :param N_points: The number of points to sample the velocity in.
    :type N_points: int
    """
    if verbose:
        print 'Sampling a VTK file: {}'.format(file)

    solver=solver or postprocess.determineSolver(file)
    if type(file) == str:
        vtkFile = vtkextract.vtkFile()
        vtkFile.solver = solver
        vtkFile.readFile(file, verbose=False)
    else:
        vtkFile = file
    
    centers = vtkFile.getCenters()
    vtkFile.getVolumes()
    bounds = vtkFile.bounds
        
    p1 = p1 or [x_c, 0, -1]
    p2 = p2 or [x_c, 0, np.max(bounds[:, 4:])]
    
    if len(p1) == 2: p1.insert(1, 0)  #if 2D data is provided
    if len(p2) == 2: p2.insert(1, 0)

    line = vtk.vtkLineSource()
    line.SetResolution(N_points - 1)
    line.SetPoint1(p1)
    line.SetPoint2(p2)
    
    probe = vtk.vtkProbeFilter()
    probe.SetInputConnection(line.GetOutputPort())
    probe.SetSourceConnection(vtkFile.outputPort)
    
    probe.Update()
    probe = probe.GetOutput()
    
    vel=probe.GetPointData().GetArray(vtkFile.velArray[solver])
    vel=np.array([vel.GetTuple(j) for j in range(probe.GetNumberOfPoints())])
    line=np.array([probe.GetPoint(j) for j in range(probe.GetNumberOfPoints())])
    
    return line, vel

def enSightTimeSteps(file):
    ts = 0
    for l in open(file).readlines():
        if 'number of steps' in l:
            ts = int(l.split()[-1])
    return ts

def sample_ensight(file, p1=None, p2=None, resolution=100, solver='TruchasEnSight', vtkFile=None):
    vtkFile = vtkextract.vtkFile()
    vtkFile.solver = solver
    vtkFile.sim='3D'
    
    ts = enSightTimeSteps(file)
    vtkFile.readFile(file, ts=ts-1)
    
    return sample_vtk(vtkFile, p1, p2, resolution, solver)
    

def cylinderReference(x, y, R, U=1):
    r = np.sqrt(x**2 + y**2)
    theta = np.arctan2(y, x)
    V_r=U*(1-R**2/r**2)*np.cos(theta)
    V_theta=-U*(1+R**2/r**2)*np.sin(theta)
    return [V_r*np.cos(theta) - V_theta*np.sin(theta), V_r*np.sin(theta) + V_theta*np.cos(theta)]
    
def cylinder(path = './', solver=None, cs='cfr'):
    if len(cs) > 1:
        for cs0 in cs:
            cylinder(path, solver, cs0)
        return
    path = os.path.abspath(path)
    solver = solver or postprocess.determineSolver(path)
    if solver == 'TruchasEnSight': solver = 'Truchas'
    U=1.
    R=0.05
    if solver == 'Gerris':
        gfsFile=glob.glob(os.path.join(path, '*.gfs'))[0]
        for l in open(gfsFile).readlines():
            if 'SurfaceBc U' in l:
                U = -float(l.split()[-1])
            if 'SolidMoving' in l:
                for i, w in enumerate(l.split()):
                    if w == 'scale':
                        R = 0.25 * float(l.split()[i+2])
            if 'Define R' in l:
                R = float(l.split()[-1])
    elif solver == 'Truchas':
        inpFile=glob.glob(os.path.join(path, '*.inp'))[0]
        for l in open(inpFile).readlines():
            l = l.split()
            if len(l) ==3 and l[0] == 'R' and l[1] == '=':
                R = float(l[2])
    elif solver == 'OpenFOAM':
        inpFile=glob.glob(os.path.join(path, 'system/snappyHexMeshDict'))[0]
        for l in open(inpFile).readlines():
            l = (l.replace(';', '')).split()
            if len(l) == 2 and l[0] == 'radius':
                R = float(l[1])
    position = { 'c':0, 'f' : -2*R, 'r' : 2*R }
    y0 = { 'c':R, 'f' : 0, 'r' : 0 }
    cs_name={ 'c':'center', 'f' : 'front', 'r' : 'rear'}
    
    fig = plt.figure()
    fig.subplots_adjust(0.15, 0.15, 0.95, None)
    if 'static' in path.lower(): 
        state = 'static'
    else:
        state = 'moving'
    if 'slip' in path:
        BC = 'slip'
    else:
        BC = 'outflow'

    plt.title('Cylinder - velocity profiles (' + ', '.join([solver, state, BC, cs_name[cs]]) + ')')
    
    r = np.linspace(y0[cs], 0.5, 100)

    ls=['-.', '--', '-', '-', '-']
    for i, caseDir in enumerate(sorted(glob.glob(os.path.join(path, '*[1-9]')))[:5]):
        vtks = postprocess.globFiles(caseDir, ext='.vtk', exit=0)
        if vtks:
            print 'Profiles:', caseDir, '*.vtk'
            line, u = sample_vtk(vtks[-1], [position[cs], y0[cs]], [position[cs], 0.5])
            if 'moving' in path:
                u = u+U     #we want relative speed
            plt.plot(u[:, 0], line[:, 2], label='Resolution '+caseDir[-1], ls=ls[i])
    
        files = glob.glob(os.path.join(caseDir, '*.ensight.CASE'))
        if files:
            file = files[0]
            print 'Profiles:', file
            ts = enSightTimeSteps(file)
            if len(glob.glob(os.path.join(caseDir, '*.ensight/P.*'))) == ts:  #check if computations are finished
                line, u = sample_ensight(file, [position[cs], y0[cs]], [position[cs], 0.5])
                if 'moving' in path:
                    u = u+U     #we want relative speed
                plt.plot(u[:, 0], line[:, 2], label='Resolution '+caseDir[-1], ls=ls[i])

    Vx, Vy = cylinderReference(position[cs], r, R)
    #plt.plot(U*(1+R**2/r**2), r, label='Reference', lw=2)
    plt.plot(Vx, r, label='Reference', c='k', lw=2.0)
    
    legend = { 'c':'upper right', 'f' : 'upper left', 'r' : 'upper left' }
    plt.legend(loc=legend[cs])
    #plt.ylim([0, line[:, 2].max()])
    #plt.ylim([0, 0.15])
    if cs == 'c':plt.xlim([1, 2])
    
    plt.ylabel('z $[m]$')
    plt.xlabel('Velocity $u_x$ $[\\frac{m}{s}]$')
    
    #name = os.path.split(path)[1]  #ex. cylinder-static, cylinder-moving
    plt.savefig('-'.join(['cylinder', solver, 'profiles', state, BC, cs]))
    plt.savefig('-'.join(['cylinder', solver, 'profiles', state, BC, cs])+'.pdf')

def get_reference_fenton(a, t0=0, x_c=15):
    if type(t0) == list:
        return np.array([r for r in imap(getReference, repeat(a), t0, repeat(x_c))])
    
    #print 'Getting a reference solution'
    subprocess.call('solitonVelocitySection '+str(a)+' 15 '+str(t0) +' '+str(x_c), shell=True)
    ref = np.loadtxt('u_section.dat')
    os.remove('u_section.dat')
    raise NotImplementedError('Reshape the result in accordance to get_reference_bim')
    return ref

def get_reference_bim(a, t0=0, x_c=0, x0=15, verbose=True):
    if type(t0) == list:
        return np.array([r for r in imap(getReferenceBIM, repeat(a), t0, repeat(x_c))])
    
    if verbose: 
        print 'Getting a reference solution for a={} from BIM data'.format(a)
   
    numRefDir = os.path.join(os.environ['HOME'], 'work/soliton/fullPotentialSolution')
    if not(os.path.exists(numRefDir)):
        sys.exit('Numerical reference directory does not exist: '+numRefDir)

    x_c = x_c - solitonVelBIM[a]*t0 - x0
    N=200
    line = (np.ones(N)*x_c, np.linspace(-1, a, N))
    
    u, ext = postprocess.readGphov(os.path.join(numRefDir, str(a), 'u'))
    v, ext = postprocess.readGphov(os.path.join(numRefDir, str(a), 'v'))
    grid_x, grid_y = np.mgrid[ext[0]:ext[1]:u.shape[1]*1j, ext[2]:ext[3]:u.shape[0]*1j]
    
    u = u.transpose()
    v = v.transpose()
    
    ux_sampled = griddata((grid_x.flatten(), grid_y.flatten()), u.flatten(), line, method='linear', fill_value=0)
    uy_sampled = griddata((grid_x.flatten(), grid_y.flatten()), v.flatten(), line, method='linear', fill_value=0)

    return np.array(line).transpose(), np.array([ux_sampled, uy_sampled]).transpose()

get_reference = get_reference_bim

def plotCrossSections(path, time=None, a=None, solver=None, x_c=25, bw=0, suffix='', difference=False):
    """
    Plots the velocity profiles along vertical crossections of a solitary wave, for a single time step.
    In previous versions multiple timesteps were supported, but is not any more.

    :param difference: Take a difference between the reference solution and numerical results
    :type difference: bool
    """

    solver = solver or postprocess.determineSolver(path)
    a = a or postprocess.determineAmplitude(path)
    
    path=os.path.abspath(path)
    dirs = glob.glob(path + '/case[1-4]/*postprocessed.dat') + glob.glob(path + '/*-[1-4]/*postprocessed.dat')
    dirs = [os.path.split(dir)[0] for dir in dirs]
    dirs.sort()
    if not (dirs): 
        print 'Empty directory:', path
        return 1
    else: 
        print dirs
    
    #correct time:
    dt = postprocess.findDt(path)
    time_step = int(round(time/dt))
    time = time_step*dt

    ref_line, u_ref = get_reference(a, time, x_c)
    y_ref = ref_line[:, 1]
    ux_ref = u_ref[:, 0]

    ux_num_list, sample_line_z_list = [], []
    
    for path in dirs:
        vtk_file = postprocess.globFiles(path, ext='.vtk', verbose=False)[time_step]
        if difference:
            p1 = [x_c, np.min(y_ref)]
            p2 = [x_c, np.max(y_ref)]
            line, vel = sample_vtk(vtk_file, p1=p1, p2=p2, N_points=len(y_ref))
            vel[:, 0] = vel[:, 0] - ux_ref
        else:
            line, vel = sample_vtk(vtk_file, x_c=x_c)
        ux_num_list.append(vel[:, 0])  #x component
        sample_line_z_list.append(line[:, 2]) #vertical component
    
    print "Plotting..."
    
    gravity = 'sourceGrav' if 'sourceGrav' in path else ''
    results = ['VelCross', "A0{0:g}".format(a*10), solver, gravity, suffix]
    results = [p for p in results if p]  #remove empty strings
    results = '_'.join(results)
    
    xlim = {0.1: [-0.03, 0.45], 0.3: [-0.3, 1.3]}
    plt.clf()
    fig = plt.figure()
    fig.subplots_adjust(0.15, 0.15, 0.95, None)
    #plt.title('Propagation of a {2:.1f} soliton, t={0:.2f} s, $\Delta x$={1:.2f}m, {3}'.format(t, x_c-15, a, solver))
    
    ls = ['-.', '--', '-']
    for j in range(len(ux_num_list)):
        case = dirs[j][-1:]
        u = ux_num_list[j]
        sample_line_z = sample_line_z_list[j]
        if bw:
            plt.plot(u, sample_line_z, label='Resolution ' + case, ls=ls[j], color='k')
        else:
            plt.plot(u, sample_line_z, label='Resolution ' + case, ls=ls[j])
    if not difference:
        plt.plot(ux_ref, y_ref, label='Reference', color='k', ls='-', lw=2)
    plt.legend(loc='best')
    
    plt.ylabel('z $[m]$')
    plt.xlabel('Velocity $[\\frac{m}{s}]$')
    plt.ylim(-1, np.max(sample_line_z))
    plt.xlim(xlim[a])
    if bw:
        results = results+'_bw'
    plt.savefig(results+'.png'.format(time))
    plt.savefig(results+'.pdf'.format(time))

    #limits for multiple timesteps:
    #v_lim1 = np.min([np.min(ux_num_list), np.min(ux_ref)])
    #v_lim2 = np.max([np.max(ux_num_list), np.max(ux_ref)])
    #plt.xlim(v_lim1, v_lim2)
    #if not(os.path.exists(results)):
    #    os.mkdir(results)
    #plt.savefig(os.path.join(results, 'velocity_{0:.2f}.png'.format(t)))
    #plt.savefig(os.path.join(results, 'velocity_{0:.2f}.pdf'.format(t)))

def plot_cross_sections_for_paper():
    x0 = 15
    time = 10.0; 
    a = 0.3
    x_c = time*solitonVel[a] + x0

    for bw in [1, 0]:
        plotCrossSections('/work/comp/foam/propagation-A03', time=time, x_c=x_c, bw=bw, suffix='calpha10')
        plotCrossSections('/work/comp/foam/propagation-A03-calpha05', time=time, x_c=x_c, bw=bw, suffix='calpha05')
        plotCrossSections('/work/comp/foam/propagation-A03-calpha20', time=time, x_c=x_c, bw=bw, suffix='calpha20')

        plotCrossSections('/work/comp/foam/propagation-A03', time=time, x_c=x_c, bw=bw)
        plotCrossSections('/work/comp/gerris/propagation-A03-sourceGrav', time=time, x_c=x_c, bw=bw)
        plotCrossSections('/work/comp/gerris/propagation-A03', time=time, x_c=x_c, bw=bw)
        plotCrossSections('/work/comp/thetis/propagation-A03', time=time, x_c=x_c, bw=bw)
        plotCrossSections('/work/comp/truchas/propagation-A03', time=time, x_c=x_c, bw=bw)

def plot_cross_sections_for_review_3():
    x0 = 15
    time = 10.0; 
    a = 0.3
    x_cross_middle = time*solitonVel[a] + x0
    delta = 1

    for x_c, suffix in [(x_cross_middle - delta, 'back'), (x_cross_middle, 'center'), (x_cross_middle + delta, 'front')]:
        for bw in [1, 0]:
            plotCrossSections('/work/comp/foam/propagation-A03', time=time, x_c=x_c, bw=bw, difference=True, suffix=suffix)
            plotCrossSections('/work/comp/gerris/propagation-A03-sourceGrav', time=time, x_c=x_c, bw=bw, difference=True, suffix=suffix)
            plotCrossSections('/work/comp/gerris/propagation-A03', time=time, x_c=x_c, bw=bw, difference=True, suffix=suffix)
            plotCrossSections('/work/comp/thetis/propagation-A03', time=time, x_c=x_c, bw=bw, difference=True, suffix=suffix)
            plotCrossSections('/work/comp/truchas/propagation-A03', time=time, x_c=x_c, bw=bw, difference=True, suffix=suffix)


if __name__ == '__main__':
    plot_cross_sections_for_review_3()
    #plot_cross_sections_for_paper()

