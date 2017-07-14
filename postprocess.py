"""The script postprocesses *.vtk files"""

import sys, os, glob, subprocess, string, re, time, shutil, ipdb
import numpy as np
import matplotlib
from ipdb import set_trace as debug

matplotlib.use('Agg', warn=False)  #it's needed for TeX in plots, but disable interactive mode in ipython, because this file is automatically imported
from matplotlib import rc
import matplotlib.pyplot as plt
import vtkextract
import numbers
import scipy.interpolate
import scipy.optimize

rc('text', usetex=True)

g = 9.80665
Ek_ref = {0.0: 1, 0.1: 0.0259711*1000*g, 0.3:0.1480683*1000*g}
Ep_ref = {0.0: 1, 0.1: 0.0249971*1000*g, 0.3:0.1343140*1000*g}
mass_ref = {0.0: 1, 0.1:0.7545816*1000, 0.3:1.3707443*1000}
h_ref = {10:{0.0: 1, 0.1:0.369, 0.3:1.275}, 30:{0.3:0.8691}}  #dicionary: h_ref[slope][a]

#used for directory issues
solvers = ['gerris', 'truchas', 'thetis', 'foam']
testCases = ['propagation-A01', 'propagation-A03', 'runup-A01', 'runup-A03', 'stillwater']
inputExt = {  'Gerris': '.gfs', 'OpenFOAM': '.foam', 'Thetis': '.don', 'Truchas': '.inp', 'TruchasEnSight': '.inp'}

def removeDuplicates(seq):
    """
    Remove duplicates from a list or a sequence
    """
    seen = set()
    seen_add = seen.add
    return [ x for x in seq if x not in seen and not seen_add(x)]

def natural_sort(l):
    """
    Required to sort the filenames; taken from StackOverflow
    """
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    return sorted(l, key = alphanum_key)
    
def readGphov(path, extrapolate=True, rescale=True, enlarge = True, N_extrapolate = 5):
    """
    Reads gphov type files (velocities and pressure from BIM/Tanaka solver)
    """
    f=open(path)
    par=f.readline().split()
    Nx = int(par[1])  #number of nodes, not cells. Nc = N-1
    Ny = int(par[2])
    dx = float(par[3])
    dy = float(par[4])
    x0 = float(par[5])
    y0 = float(par[6])
    
    u = np.loadtxt(f).flatten().reshape((Ny, Nx))
    if enlarge: #enlargement of the domain/frame to avoid picking values from outside
        add_top = 5
        add_bottom = 4
        u = np.vstack((u, -666*np.ones((add_top, Nx))))  #top
        u = np.vstack((-555*np.ones((add_bottom, Nx)), u))  #bottom
        y0 = y0-add_bottom*dy
        Ny = Ny + add_top + add_bottom
    for col in u.transpose():  #could be probably made faster, but at least it works
        mask = col<-600
        level_set = np.arange(mask.size)-np.invert(mask).sum() + 1 #subsequent positive numbers instead of -666, 0 marks the last water filled cell
        band = np.logical_and(mask, level_set<=N_extrapolate)  #a band of N_extrapolate cells, initially filled with air
        #debug()
        col[band]=col[col>-500][-1]   #all the 'air' values in a band above water surface are replaced by the last 'water' value
        col[col<-600]=0  #all the remaining velocities in the air are set to 0 (normally their value is -666)
        col[col<-500]=col[col>-500][0]   #all the values extended below the bottom are replaced by the bottom 'water' value
    extent = [x0, x0 + (Nx-1)*dx, y0, y0 + (Ny-1)*dy]  #dx * (no. of cells)
    if rescale:
        u = u*np.sqrt(g)
    f.close()
    
    return (u, extent)

def determineSolver(path='./', verbose=False):
    """
    Determines the solver used for simulation, on the basis of the path
    """
    path = os.path.abspath(path)
    if 'gerris' in path.lower(): solver = 'Gerris'
    elif 'truchas' in path.lower(): solver = 'Truchas'
    elif 'foam' in path.lower(): solver = 'OpenFOAM'
    elif 'thetis' in path.lower(): solver = 'Thetis'
    else: 
        sys.exit('Cannot determine solver by path: ' + path)
    if verbose:
        print "Solver determined by path:", solver
    
    if (os.path.splitext(path)[1] == '.inp' and 'moving_solid' in open(path).read()) or (solver=='Truchas' and 'cylinder' in path):
        solver = 'TruchasEnSight'   #special case for old version of Truchas
    return solver

def globFiles(path='./', ext='.vtk', verbose=False, exit=1):
    """
    Return a list of files with given extension. 
    
    Path does not have to be a directory. Any given path will be extended up to * + ext.
    Useful for multiple basenames in one directory - user can provide a path ending with a basename.
    """
    
    if os.path.isdir(path):
        pattern = os.path.abspath(os.path.join(path, '*' + ext))
    else:
        pattern = os.path.abspath(path + '*' + ext)
    files=glob.glob(pattern)
    
    if len(files) ==0:
        if exit:
            sys.exit("No " + pattern + " files were found")
        else:
            if verbose:
                print "No " + pattern + " files were found"
            return []
    files=natural_sort(files)
    
    #workaround for Gerris producing a vtk file in TWO first timesteps (bug)
    solver = determineSolver(files[0])
    if solver == 'Gerris' and '.0000.vtk' in files[0] and '.0001.vtk' in files[1]:
        files.pop(1)
    
    if verbose:
        print 'Found', len(files), pattern, 'files'
        for file in files:
            print os.path.split(file)[1], 
        print
    return files
    
def determineAmplitude(path='./'):  
    """
    Determines the amplitude of the solitary wave used for simulation. 
    Returns a float: 0.0, 0.1 or 0.3
    """
    path = os.path.abspath(path)
    if 'A01' in path: return 0.1
    elif 'A03' in path: return 0.3
    elif 'A00' in path: return 0.0
    else: 
        print 'Cannot determine amplitude by path: ' + path + ', return 0.0'
        return 0.0

def determineSlope(path='./'): 
    """
    Determines the slope of the beach used in runup simulation.
    Return an integer: 10 or 30 (degrees)
    """
    path = os.path.abspath(path)
    if 'slope30' in path: 
        return 30
    else: 
        return 10

def determineTestCase(path='./', verbose=False): 
    """
    Returns the name of the test case: 'Propagation ' or 'Runup ' + a
    """
    path = os.path.abspath(path)
    a = determineAmplitude(path)
    if 'propagation' in path.lower(): test = 'Propagation ' + str(a)
    elif 'runup' in path.lower(): test = 'Run up ' + str(a)
    else: 
        print 'Cannot determine test case by path: ' + path
        test = 'Stillwater'
    if verbose:
        print "Test case determined by path:", test
    return test

def caseDetailsFromPath(path='./', verbose=False): 
    """
    Collects the information given by the directory structure.
    Returns a string.
    """
    path = os.path.abspath(path)
    path = path.split('/')
    for dir in ['', 'work', 'comp'] + solvers:
        if dir in path:
            path.remove(dir)
    details = '/'.join(path)
    if verbose and path:
        print "Test case determined by path:", details
    return details

def caseDetailsFromAppendix(path='./', verbose=False): 
    """
    Collects the information from the name of the test directory.
    Returns a string.
    """
    path = os.path.abspath(path)
    path = os.path.split(path)[-1]
    details=''
    if len(path.split('-A0'))==2:
        details = path.split('-')[2:]
        details = '_'.join(details)
    if verbose:
        print "Test case determined by path:", details
    return details

def findBasename(file, solver=None, verbose=False):
    """
    Finds a basename of the test case, given an input file
    """
    file = os.path.abspath(file)
    solver = solver or determineSolver(file)
    filename = os.path.split(file)[1]
    if solver in ['Truchas', 'Gerris']:
        basename = filename.split('.')[0]
    elif solver in ['Thetis', 'OpenFOAM']:
        basename = os.path.splitext(filename)[0].split('_')[0]
    elif solver == 'TruchasEnSight':
        basename = os.path.splitext(os.path.splitext(filename)[0])[0]  #.ensight.CASE file
    
    if verbose:
        print 'Basename:', basename
    return basename
    
def findDt(path, basename=None, solver=None):
    """
    Determines a timestep used in a simulation
    """
    path = os.path.abspath(path)
    solver = solver or determineSolver(path)
    ext = inputExt[solver] 
    if basename:
        filename = os.path.join(path, basename + ext)
    else:
        filename = glob.glob(os.path.join(path, '*' + ext))[0]

    if solver == 'Gerris':
        lines = open(filename).readlines()
        lines = [re.split('[ ={}]+', l) for l in lines]
        lines = [l for l in lines if 'GfsOutputSimulation' in l and 'step' in l and 'VTK' in l]
        if lines:
            l = lines[0]
            i = l.index('step')
            dt = float(l[i + 1])
        else:
            sys.exit('Cannot find a timestep information in a Gerris .gfs file.')
    elif solver in ['Truchas', 'TruchasEnSight']:
        lines = open(filename).readlines()
        lines = [re.split('[ ={}, ]+', l) for l in lines]
        lines = [l for l in lines if 'Output_Dt' in l]
        if lines:
            l = lines[0]
            dt = float(l[2])
        else:
            sys.exit('Cannot find a timestep information in a Truchas .inp file.')
    elif solver == 'Thetis':
        lines = open(filename).readlines()
        lines = [re.split('[ ]+', l) for l in lines]
        l = [l for l in lines if l[0] == 'PAS_DE_TEMPS' and l[1] == 'NAVIER'][0]
        dt = float(l[2].replace('D0', ''))
        l = [l for l in lines if l[0] == 'IMPRESSION' and l[1] == 'FREQUENCE'][0]
        dt = dt*int(l[2])
    elif solver == 'OpenFOAM':
        lines = open(os.path.join(path, 'system/controlDict')).readlines()
        lines = [re.split('[ ;]+', l) for l in lines]
        l = [l for l in lines if l[0] == 'writeInterval'][0]
        dt = float(l[1])
    else:
        sys.exit('Wrong solver (findDt function)')
    return dt


def extractTestData(data='./', outPath='./', label=None, f='output_data.txt', verbose=True):
    """
    Extracts essential data from *_postprocessed.dat files. 
    The 'data' parameter is expected to be either a list of paths to *_postprocessed.dat files, 
    or a single string giving path to the main test directory.
    Saves the extracted data into the f file, saved in outPath.
    """
    if type(data) == str:
        data = glob.glob(os.path.join(data, '*[1-9]/*postprocessed.dat'))
    if not data:
        print "extractTestData: No data files provided, or no '*/*postprocessed.dat' files found"
        return 1
    data = [os.path.abspath(d) for d in data]
    data = sorted(data)
    dirs = [os.path.split(d)[0] for d in data]
    label = label or os.path.split(dirs[0])[1][:-2]
    prefix=os.path.commonprefix(dirs)
    if verbose:
        print 'Extracting output data:', prefix, [d[len(prefix):] for d in dirs]
    f=open(os.path.join(outPath, f), 'w')
    a = determineAmplitude(dirs[0])
    slope = determineSlope(dirs[0])
    P = [] #container for all postprocessed data
    for d in data:
        P.append(np.loadtxt(d))
    f.write('#' + label + ': ' + prefix + ' ' + str([d[len(prefix):] for d in dirs]) + '\n')
    for dir in dirs:
        t = open(glob.glob(os.path.join(dir, '*timing.dat'))[0]).readline().split()[0]
        f.write('{0:.0f}\t'.format(np.abs(float(t))))
    f.write('#Timing\n')
    
    for p in P:
        f.write('{0:.1f}\t'.format(p[-1, 1]))
    f.write('#final kinetic energies\n')
    for p in P:
        f.write('{0:.3f}\t'.format((p[-1, 1])/Ek_ref[a]))
    f.write('#final kinetic energies normalized with {0}\n'.format(Ek_ref[a]))
    
    if 'runup' in dirs[0]:
        for p in P:
            f.write('{0:.3f}\t'.format(np.max(p[p[:, 0]<=8.][:, -1])))
        f.write('#maximum height (up to 8th second)\n')
        for p in P:
            f.write('{0:.3f}\t'.format(np.max(p[p[:, 0]<=8.][:, -1])/h_ref[slope][a]))
        f.write('#maximum height normalized with h_ref={0}\n'.format(h_ref[slope][a]))
        for p in P:
            f.write('{0:.3f}\t'.format(np.max(p[p[:, 0]<=8.][:, -2])))
        f.write('#maximum height (up to 8th second) - single fluid region\n')
        for p in P:
            f.write('{0:.3f}\t'.format(np.max(p[p[:, 0]<=8.][:, -2])/h_ref[slope][a]))
        f.write('#maximum height normalized with h_ref={0} - single fluid region\n'.format(h_ref[slope][a]))
    else:
        for p in P:
            f.write('{0:.3f}\t'.format(p[-1, -1]))
        f.write('#final height\n')
        for p in P: 
            f.write('{0:.3f}\t'.format(p[-1, -1]/a))
        f.write('#final height normalized\n')
        for d in data:
            x = np.loadtxt(glob.glob(os.path.join(os.path.split(d)[0], '*center_max_height.dat'))[0])
            f.write('{0:.4f}\t'.format(x[-1]-15))
        f.write('#distance travelled by the wave; final position computed from maximum height\n')
        for d in data:
            x = np.loadtxt(glob.glob(os.path.join(os.path.split(d)[0], '*center_full_width.dat'))[0])
            f.write('{0:.4f}\t'.format(x[-1]-15))
        f.write('#distance travelled by the wave; final position computed from the center of full width at half maximum\n')
        
    f.close()
    return 0
    
class plotMaker(object):
    """
    A general object used for plotting variables
    """
    
    def __init__(self, t, file_name, xlabel='Time $[s]$', out_dir='./', loc='best', text_pos=[0, 0.01], info='', a=0.1, slope=10, ext=['png'], title=""):
        #t should be a list of numpy arrays, or a numpy array itself
        self.t = np.array(t)
        if self.t.ndim == 1:  #in case a single set of data is provided
            self.t=[self.t, ]
        self.file_name = file_name + '.{ext}'
        self.xlabel = xlabel
        self.title = title
        self.info = info
        self.loc = loc
        self.text_pos = text_pos
        self.out_dir = out_dir
        self.slope = slope
        self.a = a
        self.ext = ext if type(ext)==list else [ext, ]

    def plotReference(self):
        """
        Adds a reference plot to the current plot
        """
        if not self.a:
            ref = np.array([[0, 0], [100, 0]])
        else:
            ref = np.loadtxt(os.path.join(os.getenv('HOME'), 'work/soliton/runup_reference/convergence/runup_max_plot_{:.1f}_slope{}.txt'.format(self.a, self.slope)))
        plt.plot(ref[:, 0], ref[:, 1], label='Reference', c='k', lw=2)

    def plotEvolution(self, var, name, name_long='new plot', figlabel='Resolution', ylabel='', loc='', reference=False):
        """
        Plots evolution in time of the given variable.
        """
        plt.clf()
        var = np.array(var)
        if var.ndim == 1:  #in case a single set of data is provided
            var=[var, ]
        ls=['-.', '--', '-']
        kwargs = {}
        if 'bw' in self.file_name:
            kwargs['color'] = 'k'
        for i in range(len(self.t)):
            plt.plot(self.t[i], var[i], label='{0} {1}'.format(figlabel, i + 1), ls=ls[i], **kwargs)
        plt.title(self.title.format(name_long))
        plt.xlabel(self.xlabel)
        plt.ylabel(ylabel)

        if self.info:
            plt.figtext(self.text_pos[0], self.text_pos[1], self.info)
        plt.xlim(np.round(np.min(self.t[0])), np.round(np.max(self.t[0])))
        if reference: self.plotReference()
        plt.legend(loc=loc or self.loc)
        for e in self.ext:
            plt.savefig(os.path.join(self.out_dir, self.file_name.format(name=name, ext=e)))

    def plotVar_backup(self, var, fileName, label='Resolution', title='Plot', ylabel='', loc='', text=[0, 0, ''], reference=False, slope=10):
        plt.clf()
        if (type(var)==list and isinstance(var[0], (float, np.float32, np.float64))) or (type(var)==np.ndarray and var.ndim == 1):  #in case a single set of data is provided
            var=[var, ]
            t=[t, ]
        ls=['-.', '--', '-']
        for i, var_i in enumerate(var):
            if 'bw' in fileName:
                plt.plot(self.t[i], var_i, label='{0} {1}'.format(label, i + 1), ls=ls[i], c='k')
            else:
                plt.plot(self.t[i], var_i, label='{0} {1}'.format(label, i + 1), ls=ls[i])
        #plt.title(title)
        plt.xlabel(self.xlabel)
        plt.ylabel(ylabel)
        #plt.figtext(*self.text)
        plt.xlim(np.round(np.min(t[0])), np.round(np.max(t[0])))
        if reference: plotReference(plt, reference, slope)
        plt.legend(loc=loc or self.loc)
        plt.savefig(os.path.join(out_dir, fileName + '.png'))
        plt.savefig(os.path.join(out_dir, fileName + '.pdf'))



def plotCollections(test_dir='./', test=None, solver=None, plots=[], out_dir = '', bw=0, verbose=True, ending='', remove_old=False):
    """
    Plots collections of data from the given test directory. 
    """
     
    #import matplotlib.pyplot as plt
    test_dir = os.path.abspath(test_dir)
    out_dir = out_dir or test_dir
    data = glob.glob(os.path.join(test_dir, '*[1-5]', '*postprocessed.dat'))
    data = sorted(data)
    dirs = [os.path.split(data0)[0] for data0 in data]
    if data:
        test = test or determineTestCase(dirs[0])
        runup = True if 'run up' in test.lower() else False
        if not plots:
            if runup:
                plots = ['height', 'mass']
            else:
                plots = ['Ek', 'Ep', 'Et', 'height', 'mass']
        extractTestData(data, test_dir, verbose=False)
        if verbose:
            prefix=os.path.commonprefix(dirs)
            print 'Plotting collections:', prefix, [d[len(prefix):] for d in dirs]
        solver = solver or determineSolver(dirs[0])
        a = determineAmplitude(dirs[0])
        slope = determineSlope(dirs[0])
        
        plotLabel='Resolution '
        
        t, Ek, Ep, Et, Vmax, mass, z, z_largest=[[] for i in range(8)]
        for d in data:
            f=np.loadtxt(d)
            if runup:  #xlim of the plot
                N = np.sum(f[:, 0]<=11)
                f = f[:N]
            #['t', 'Ek', 'Ep', 'mass', 'massWave', 'maxVel', 'Xmax', 'Zmax_largest', 'Zmax'])
            #   0   1     2   3         4         5        6       7            8
            t.append(f[:, 0])
            Ek.append(f[:, 1] /Ek_ref[a] )
            Ep.append(f[:, 2] /Ep_ref[a] )
            Et.append((f[:, 1] + f[:, 2]) /(Ep_ref[a] + Ek_ref[a]) )
            Vmax.append(f[:, -4])
            mass.append(f[:, -5] /mass_ref[a] )
            z.append(f[:, -1])
            z_largest.append(f[:, -2])
        
        rotation = 'rotated' if 'rotated' in dirs[0] else ''
        cart = 'stretched' if 'stretched' in dirs[0] else ''
        droplets = 'remove' if 'remove' in dirs[0] else ''
        snapped = 'snapped' if 'snapped' in dirs[0] else ''
        layers = 'layers' if 'layers' in dirs[0] else ''
        gravity = 'sourceGrav' if 'sourceGrav' in dirs[0] else ''
        slopeLabel = 'slope30' if slope==30 else ''
        bw = 'bw' if bw else ''

        #file_name = [test.replace('Run up', 'Runup').replace(' ', '_').replace('.', ''), slopeLabel, label, solver, rotation, snapped, layers, cart, gravity, droplets, ending, bw]
        file_name = [test.replace('Run up', 'Runup').replace(' ', '_').replace('.', ''), '{name}', solver, caseDetailsFromAppendix(test_dir), bw]
        file_name = [b for b in file_name if b]  #remove empty
        file_name = '_'.join(file_name)

        info = '' #caseDetailsFromPath(os.path.split(dirs[0])[0]).replace('_', '\_')
        title ='' #'' "{test} - {{}} ({solver})".format(test = test, solver = solver)
        if remove_old:
            subprocess.call('rm {0}/*.png {0}/*.pdf -f'.format(out_dir), shell=True)
        plot_maker=plotMaker(t, file_name, ext=['png', 'pdf'], xlabel='Time $[s]$', info=info, out_dir=out_dir, a=a, slope=slope, title=title)
        if 'Ek' in plots:       plot_maker.plotEvolution(Ek, 'EK', name_long='kinetic energy', ylabel='Kinetic energy $[J]$')
        if 'Ep' in plots:       plot_maker.plotEvolution(Ep, 'EP', name_long='potential energy', ylabel='Potential energy $[J]$')
        if 'Et' in plots:       plot_maker.plotEvolution(Et, 'ET', name_long='total energy', ylabel='Total energy $[J]$')
        if 'mass' in plots:     plot_maker.plotEvolution(mass, 'Mass', name_long='mass', ylabel='Mass $[kg]$')
        if 'Vmax' in plots:     plot_maker.plotEvolution(Vmax, 'MaxVel', name_long='maximum velocity', ylabel='Velocity $[\\frac{m}{s}]$')
        if 'height-multi' in plots:   plot_maker.plotEvolution(z, 'MultiHeight', name_long='maximum height in the whole domain', ylabel='Height $[m]$', reference=runup)
        if 'height' in plots:plot_maker.plotEvolution(z_largest, 'Height', name_long='maximum height', ylabel='Height $[m]$', reference=runup)

    else:
        print 'plotCollections: No output directories found in', test_dir
        return 1


def plotErrors(dirs, labels, basename, test='test', out_dir = './'):
    #Function expects a list of directories with output_data.txt files
    errors=[]
    plt.clf()
    a = determineAmplitude(dirs[0])
    for dir in dirs:
        print dir
        if extractTestData(dir, verbose=False):
            errors.append([0, 0, 0])
            continue
        print dir
        f=np.loadtxt(os.path.join(dir, 'output_data.txt'))
        if 'runup' in dir:
            errors.append((1 - f[6, :])*100)
            plt.ylabel('$\\varepsilon_N [\\%]$')
            #plt.title('Run up {0:.1f} height error - {1} domain'.format(a, 'rotated' if 'rotated' in dir else 'horizontal'))
            filename='RunupHeight.png'
        else: #propagation
            errors.append(f[2, :])
            plt.ylabel('Energy decrease $[\%]$')
            #plt.title(test + ' - Energy decrease')
            filename='EnergyDecrease.png'
    plotLabel=''
    errors = np.array(errors)
    O = errors[:, 0].max()
    markers = ['s', 'o', 'D', '^']
    for i in range(len(dirs)):
        if 'bw' in basename:
            plt.plot(np.arange(1, 4), errors[i], label=labels[i], marker=markers[i], ms=10, c='k')
        else:
            plt.plot(np.arange(1, 4), errors[i], label=labels[i], marker=markers[i], ms=10)
    
    if 'bw' in basename:
        plt.plot(np.arange(1, 4), [O, O/2, O/4], label='O(1)', marker='s', mfc='none', ms=10, c='k', ls='--')
    else:
        plt.plot(np.arange(1, 4), [O, O/2, O/4], label='O(1)', marker='s', mfc='none', ms=10, c='grey', ls='--')

    plt.xlabel('Resolution')
    plt.xlim((0.5, 3.5))
    plt.xticks(np.arange(1, 4))
    plt.legend(loc='best')

    plt.savefig(os.path.join(out_dir, basename + '.png'))
    plt.savefig(os.path.join(out_dir, basename + '.pdf'))
    
def compareAllErrors(out_dir = './'):
    compDir='/work/comp/'
    dirs = []; labels=[]

    dirs.append('gerris/'); labels.append('Gerris')
    dirs.append('foam/'); labels.append('OpenFOAM')
    dirs.append('thetis/'); labels.append('Thetis')
    dirs.append('truchas/'); labels.append('Truchas')

    #for test in ['propagation-A01', 'propagation-A03', 'runup-A01', 'runup-A03']:
    for test in ['runup-A01', 'runup-A03']:
        print test
        dirsFull=[os.path.join(compDir, dir, test) for dir in dirs]
        basename = test + '_errors_'
        plotErrors(dirsFull, labels, basename, test)

def compareRunupErrors(compDir='/work/comp/', out_dir = './', bw=0):
    dirs = []; labels=[]
    
    labels = ['Gerris', 'OpenFOAM', 'Thetis', 'Truchas']
    dirs = ['gerris/runup-A0{0}{1}', 'foam/runup-A0{0}{1}', 'thetis/runup-A0{0}{1}', 'truchas/runup-A0{0}{1}']

    for A in [1, 3]:
        for rotated in ['', '-rotated']:
            dirsTest = [os.path.join(compDir, d).format(A, rotated) for d in dirs]
            
            basename = 'Runup_errors-0{0}{1}'.format(A, rotated)
            if bw:
                basename = basename + '_bw'
            plotErrors(dirsTest, labels, basename)


def postprocessDirs(dirs='./', solver=None):
    """
    Helper function, running postprocess function in several directories
    dirs parameter is either a list of directories, or a path to the main test directory
    """
    if type(dirs)==str:  #in case a path is provided to the main catalogue
        dirs = os.path.abspath(dirs)
        dirs = [dir for dir in glob.glob(dirs + '/*[1-5]') if glob.glob(dir + '/*.vtk')]
        if not dirs:
            print 'postprocessDirs: No output directories found'
            return 1
    for dir in sorted(dirs):
        postprocess(dir, solver=solver)

def saveData(dir='./', basename=None, var=[], labels=[], data_1_file={}, data_sep_files={}):
    """
    Saves the data extracted from VTK files into text files (e.g. *_postprocessed.dat)
    """
    dir = os.path.abspath(dir)
    basename = basename or findBasename(dir)
    f=open(os.path.join(dir, basename + '_postprocessed.dat'), 'w')
    f.write(('#' + '{:^11s}'*len(labels) + '\n').format(*labels))
    var = np.vstack(var).transpose()
    np.savetxt(f, var, fmt='%10.6f')
    f.close()
    for key, value in data_sep_files.iteritems():
        f_name=os.path.join(dir, '{}_data_{}.dat'.format(basename, key))
        #var = np.vstack(value).transpose()
        np.savetxt(f_name, value, fmt='%10.6f')


def postprocess(path, solver=None):
    """
    Postprocesses output directory of a simulation.
    Use it for 2 phase simulations: propagation of a solitary wave, run up of a solitary wave, 
    still water test.
    Extracts data from all VTK files, computes and plots kinetic energy, maximum velocity etc.
    """
    
    path = os.path.abspath(path)
    solver = solver or determineSolver(path)
    test = determineTestCase(path)
    amp = determineAmplitude(path)
    
    files = globFiles(path, ext='.vtk', verbose=False, exit=0)
    if not files:
        return 1
    basename = findBasename(files[0], solver)
    
    out_dir = os.path.abspath(path) if os.path.isdir(path) else os.path.split(path)[0]

        
    dt = findDt(out_dir, basename, solver)

    vtkFile = vtkextract.vtkFile()
    vtkFile.solver = solver

    Ek, Ep, Vmax, max_x, max_z, max_z_largest, mass, mass_wave = [], [], [], [], [], [], [], []
    variable_names = ['t', 'center_full_width', 'full_width', 'center_max_height']
    data = {name: [] for name in variable_names}
    
    print 'postprocess: ', len(files), 'files in', out_dir, ': ', basename + '.#.vtk', 
    for i, file in enumerate(files):
        print i, 
        sys.stdout.flush()
        vtkFile.readFile(file, verbose=False)
        vel=vtkFile.getVelocity()
        u=np.sqrt(vel[:, 0]**2 + vel[:, 1]**2 + vel[:, 2]**2)

        vol=vtkFile.getVolumes()
        #ipdb.set_trace()
        water=vtkFile.getWaterFraction()
        rho=water*1000 + (1-water)*1
        Ekk = u**2*rho*vol/2
        Ek.append( np.sum(Ekk[water>0.005]) )
        Vmax.append(np.max(u*water))
        mass.append(np.sum(vol*rho))
        mass_wave.append( (np.sum(vol*water) - np.max(vtkFile.bounds[:, 1]))*1000)
        
        if vtkFile.time >=0:
            data['t'].append(vtkFile.time)
        else:
            data['t'].append(i*dt)
        
        contour = vtkFile.getContour()
        #sorting the x values (will spoil the results for complicated shapes, but works fine for propagation and run up:
        arg=np.argsort(contour[:, 0], axis=0)
        contour = contour[arg]
        contour = contour[contour[:, 2]>=-0.05]
        
        maxWaterLevel = np.max(contour[:, 2])
        centerX = np.mean(contour[contour[:, 2]==maxWaterLevel, 0])
        max_x.append(centerX)
        max_z.append(maxWaterLevel)
        data['center_max_height'].append(centerX)
        
        contour = vtkFile.getContour(largestRegion=True)
        arg=np.argsort(contour[:, 0], axis=0)
        contour = contour[arg]
        contour = contour[contour[:, 2]>=-0.05]
        max_z_largest.append(np.max(contour[:, 2]))
        
        #calculation of potential energy - integration of the interface
        dx = contour[1:, 0] - contour[:-1, 0]
        height = (contour[:-1, 2] + contour[1:, 2])/2
        Ep1=height*height*dx*1000*g /2
        Ep1 = np.sum(Ep1)
        Ep.append(Ep1)

        if 'propagation' in path.lower():

            #calculation of the full width at half maximum
            #making the contour 1D:
            indices = np.unique(contour[:, 0], return_index=True)[1]
            contour = contour[indices]
            contour = contour[:, [0, 2]]
            f = scipy.interpolate.interp1d(contour[:, 0], contour[:, 1] - float(amp)/2, bounds_error=False, fill_value=0)
            a = scipy.optimize.brentq(f, np.min(contour[:, 0]), centerX)
            b = scipy.optimize.brentq(f, centerX, np.max(contour[:, 0]))
            data['full_width'].append( b-a )
            data['center_full_width'].append( (a + b)/2 )



    print
    
    #for Gerris, the centers and bounds are moved by 1 in z direction by vtkextract.py
    saveData(out_dir, basename, var=[data['t'], Ek, Ep, mass, mass_wave, Vmax, max_x, max_z_largest, max_z], 
             labels=['t', 'Ek', 'Ep', 'mass', 'massWave', 'maxVel', 'Xmax', 'Zmax_largest', 'Zmax'], 
             data_sep_files=data)
    
    print "Plotting..."
    title = '{test} - {{}} ({solver})'.format(test=test, solver=solver)
    maker=plotMaker(data['t'], basename + ".{name}", ext=['png', 'pdf'], xlabel='Time $[s]$', a=determineAmplitude(path), out_dir=out_dir, slope=determineSlope(path), title=title)
    maker.plotEvolution(Ek, 'EK', name_long='kinetic energy', ylabel='Kinetic energy $[J]$', )
    maker.plotEvolution(Ep, 'EP', name_long='potential energy', ylabel='Potential energy $[J]$', )
    maker.plotEvolution(Vmax, 'MaxVel', name_long='maximum velocity', ylabel='Velocity $[\\frac{m}{s}]$', )
    runup = True if 'run up' in test.lower() else False
    maker.plotEvolution(max_z, 'Height-multi', name_long='maximum height', ylabel='Height $[m]$', reference=runup)
    maker.plotEvolution(max_z_largest, 'Height', name_long='maximum height of largest region', ylabel='Height $[m]$', reference=runup)


def plotChosenCollections(outDir0='', postprocess=False, remove_old=False):
    
    solverDirs=['gerris', 'foam', 'thetis', 'truchas']
    #testDirs=['propagation-A01', 'propagation-A03', 'runup-A01', 'runup-A03']
    #solverDirs=['foam', 'gerris']
    testDirs=['propagation-A01', 'propagation-A03']
    #testDirs=['runup-A01', 'runup-A03']
    
    out_dir = ''
    if outDir0:
        out_dir = os.path.abspath(outDir0)
    if remove_old:
        subprocess.call('rm ' + out_dir + '/*.png -f', shell=True)
    labels = []
    data = []
    for testDir in testDirs:
        for solverDir in solverDirs:
            dirs = glob.glob(os.path.join('/work/comp', solverDir, testDir + '*'))
            dirs = [dir for dir in dirs if 'test' not in dir]
            for dir in dirs:
                if not outDir0:
                    out_dir = dir
                if postprocess: postprocessDirs(dir)
                if 'runup' in testDir:
                    plots=['height']#, 'mass']
                else:
                    plots=['Et', 'height', 'Vmax']#, 'mass']
                #err = plotCollections(dirs=dir, out_dir = out_dir, plots=plots, verbose=True, ending='', bw=1)
                err = plotCollections(dirs=dir, out_dir = out_dir, plots=plots, verbose=True, ending='', remove_old=remove_old)
                
                
                if not err:
                    labels.append(dir.replace('/work/comp/', ''))
                    A = np.loadtxt(os.path.join(dir, 'output_data.txt'))
                    if len(A[0]) < 3:
                        A = np.hstack((A, A[:, -1:]))
                        A[0, -1] = 999999
                    data.append(A)
    
    f = open('output_data_all.txt', 'w')
    f.write('#{:18}|{:24}\n'.format('Normalized run up height', 'Timing'))
    #ordering of the results, for easy coping to the article
    zipped = zip(data, labels)
    #zipped = filter(lambda x:'runup' in x[1], zipped)  #only runup
    #zipped = filter(lambda x:not ('gerris' in x[1] and 'remove' not in x[1]), zipped)
    zipped = filter(lambda x:not ('runup' in x[1] and 'gerris' in x[1] and 'remove' in x[1]), zipped)
    zipped = filter(lambda x:'snapped' not in x[1] and 'layers' not in x[1], zipped)
    zipped = sorted(zipped, key=lambda x:'sourceGrav' in x[1])  #source term gravity in Gerris before the reduced gravity
    zipped = sorted(zipped, key=lambda x:'A03' in x[1])  #'A01' before 'A03'
    zipped = sorted(zipped, key=lambda x:'rotated' in x[1])  #'rotated' after 'horizontal'
    zipped = sorted(zipped, key=lambda x:'slope30' in x[1])  #30 degrees after 10
    zipped = sorted(zipped, key=lambda x:'propagation' in x[1])  #runup before propagation
    data, labels = zip(*zipped)
    #data = np.array(data)
    groups = np.array([4, 4, 4, 4, 4, 5, 5])
    ig = 0
    for i, d in enumerate(data):
        if not i or not (i % groups[0:ig].sum()):  #I need this to make groups among results with given time normalization
            time_ref= np.min(np.array(data[i:i + groups[ig]])[0, :], axis=0)
            f.write(('{:8.0f}\t'*3).format(*time_ref) + '    #time reference\n')
            ig = ig + 1
        f.write(('{:6.3f}\t'*3).format(*d[-1, :]) + ('{:8.3g}\t'*3).format(*(d[0, :]/time_ref)) + ('{:8.0f}\t'*3).format(*d[0, :]) + '    #{}\n'.format(labels[i]))
    f.close()



def postprocess1phase(path, solver=None):
    """
    Postprocesses output directory of a solver.
    Use it for 1 phase simulations: moving cylinder etc.
    Looks for all the timesteps, computes and plots kinetic energy, maximum velocity, 
    timesteps.
    """
    
    path = os.path.abspath(path)
    if os.path.isdir(path):
        out_dir = os.path.abspath(path)
    else:
        out_dir = os.path.split(path)[0]
    solver = solver or determineSolver(path)
        
    t, Ek, Vmax, Vx_max, P=[], [], [], [], []
    rho = 1000
        
    if solver == 'TruchasEnSight':
        file = globFiles(path, ext='.CASE')[0]
        basename = findBasename(file, solver)

        dt = findDt(out_dir, basename, solver)

        vtkFile = vtkextract.vtkFile()
        vtkFile.solver = solver
        vtkFile.sim='3D'
        
        ts = 0
        for l in open(file).readlines():
            if 'number of steps' in l:
                ts = int(l.split()[-1])
        
        print 'postprocess1phase: ', ts, ' ensight CASE files in', out_dir, ': ', 
        for i in range(ts):
            print i, 
            sys.stdout.flush()
            vtkFile.readFile(file, ts=i)
            #centers = vtkFile.getCenters()
            vel=vtkFile.getVelocity()
            u=np.sqrt(vel[:, 0]**2 + vel[:, 1]**2 + vel[:, 2]**2)
            vol=vtkFile.getVolumes()
            
            Ek.append(np.sum(u**2*rho*vol/2))
            Vmax.append(np.max(u))
            Vx_max.append(np.max(vel[:, 0]))
            p = vtkFile.getPressure()
            P.append(p.max())
            
            if vtkFile.time >=0:
                t.append(vtkFile.time)
            else:
                t.append(i*dt)
        print
    else:
        files = globFiles(path, ext='.vtk', verbose=False)
        basename = findBasename(files[0], solver)
            
        dt = findDt(out_dir, basename, solver)

        vtkFile = vtkextract.vtkFile()
        vtkFile.solver = solver

        print 'postprocess1phase: ', len(files), 'files in', out_dir, ': ', basename + '.#.vtk', 
        for i, file in enumerate(files):
            print i, 
            sys.stdout.flush()
            vtkFile.readFile(file, verbose=False)
            #centers = vtkFile.getCenters()
            vel=vtkFile.getVelocity()
            u=np.sqrt(vel[:, 0]**2 + vel[:, 1]**2 + vel[:, 2]**2)
            vol=vtkFile.getVolumes()
            
            Ek.append(np.sum(u**2*rho*vol/2))
            Vmax.append(np.max(u))
            Vx_max.append(np.max(vel[:, 0]))
            #p = vtkFile.getPressure()
            #P.append(p.max())
            
            if vtkFile.time >=0:
                t.append(vtkFile.time)
            else:
                t.append(i*dt)
        print
    
    t = np.array(t)
    Ek = np.array(Ek)
    Vmax = np.array(Vmax)
    #P = np.array(P)
    Vx_max = np.array(Vx_max)
    
    saveData(out_dir, basename, var=[t, Vmax, Ek, Vx_max], labels=['t', 'Vmax', 'Ek', 'Vx_max'])

if __name__ == '__main__':
    #plotChosenCollections(outDir0='', postprocess=False, remove_old=False)
    plotCollections()

